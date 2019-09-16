//#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <regex>
#include <map>
#include <memory>
#include <boost/process/system.hpp>
#include "../src/LDServer.h"
#include "../src/ScoreServer.h"
#include "RareMetal.h"
#include "../src/Phenotypes.h"
#include "../src/Mask.h"
#include "../src/ScoreCovarianceRunner.h"
#include "../src/Regression.h"
#include <armadillo>
#include <cereal/external/rapidjson/document.h>

using namespace std;

const uint32_t MAX_UINT32 = numeric_limits<uint32_t>::max();
const uint32_t INITIAL_RESULT_SIZE = 10000000;

class LDServerTest: public::testing::Test {
protected:
    static string hostname;
    static int port;
    static boost::process::child redis_server;
    static redisContext* redis_cache;

    virtual ~LDServerTest() {}

    static void SetUpTestCase() {
        ifstream config_file("redis-connection.txt");
        string hostname;
        string port;
        getline(config_file, hostname);
        getline(config_file, port);
        LDServerTest::port = stoi(port);
        LDServerTest::hostname = hostname;

        auto env_redis = getenv("REDIS_BIN");
        string redis_bin = "../bin/redis-server"; // default path to search
        if (env_redis != nullptr) {
            cout << "Found REDIS_BIN envvar, using it for redis-server path" << endl;
            redis_bin = static_cast<string>(env_redis);
        }
        else {
            cout << "No REDIS_BIN found, trying to use default path of ../bin/redis-server" << endl;
        }

        string command(redis_bin + " --port " + to_string(LDServerTest::port) + " --bind " + LDServerTest::hostname + " --daemonize no --save \"\"");
        LDServerTest::redis_server = boost::process::child(command);
        this_thread::sleep_for(chrono::seconds(3));
        redis_cache = redisConnect(LDServerTest::hostname.c_str(), LDServerTest::port);
    }

    static void TearDownTestCase() {
        if (redis_cache != nullptr) {
            redisFree(redis_cache);
        }
        LDServerTest::redis_server.terminate();
    }

    virtual void SetUp() {

    }

    virtual void TearDown() {
        redisReply* reply = (redisReply*)redisCommand(redis_cache, "FLUSHDB");
        assert(reply != nullptr);
        freeReplyObject(reply);
    }

    void load_region_goldstandard(const string &file, map<string, double> &values) {
        ifstream input_file(file);
        string line;
        vector<string> tokens;
        auto separator = regex("[ \t]");
        getline(input_file, line); // skip header;
        while (getline(input_file, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
            values.emplace(tokens.at(1) + "_" + tokens.at(2), stod(tokens.at(4)));
            tokens.clear();
        }
    }

    void load_variant_goldstandard(const string &file, map<string, double> &values) {
        ifstream input_file(file);
        string line;
        vector<string> tokens;
        auto separator = regex("[ \t]");
        getline(input_file, line); // skip header;
        while (getline(input_file, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
            values.emplace(tokens.at(1) + "_" + tokens.at(3), stod(tokens.at(5)));
            tokens.clear();
        }
    }

    shared_ptr<map<string, double>> load_variant_freqs(const string &file) {
        auto variant_freq = make_shared<map<string, double>>();
        ifstream input_file(file);
        string line;
        vector<string> tokens;
        auto separator = regex("[ \t]");
        string variant;
        double freq;
        uint64_t i = 0;
        while (getline(input_file, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
            variant = tokens.at(2);
            freq = stod(tokens.at(7));

            if ((freq < 0) || (freq > 1)) {
                throw std::domain_error("Allele frequency for variant " + variant + " is invalid: " + to_string(freq));
            }

            if (variant.empty() || (variant == ".") || (variant == "NA")) {
                throw std::domain_error("Invalid variant in frequency file on row " + to_string(i));
            }

            variant_freq->emplace(variant, freq);
            tokens.clear();
            i++;
        }
        return variant_freq;
    }

    auto load_raremetal_scores(const string &file) {
      return make_shared<RareMetalScores>(file);
    }

    void load_raremetal_covariance(const string &file, map<string, double> &values) {
        ifstream input_file(file);
        string line;
        vector<string> tokens;
        vector<uint64_t> positions;
        vector<double> cov;
        auto line_separator = regex("[ \t]");
        auto field_separator = regex(",");
        string comment = "#";

        getline(input_file, line); // skip header;
        while (getline(input_file, line)) {
            if (std::equal(comment.begin(), comment.end(), line.begin())) {
              continue;
            }

            // Split line and insert into tokens
            copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

            string chrom(tokens.at(0));
            string ref_pos(tokens.at(1));

            // Load the positions on this row
            transform(
              sregex_token_iterator(tokens.at(2).begin(), tokens.at(2).end(), field_separator, -1),
              sregex_token_iterator(),
              back_inserter(positions),
              [](const string &str) { return stoi(str); }
            );

            // Load the covariance values on this row
            transform(
              sregex_token_iterator(tokens.at(3).begin(), tokens.at(3).end(), field_separator, -1),
              sregex_token_iterator(),
              back_inserter(cov),
              [](const string &str) { return stod(str); }
            );

            // Store values to map
            string pair_key;
            double pair_cov;
            for (int index = 0; index < positions.size(); index++) {
              pair_key = ref_pos + "_" + to_string(positions[index]);
              pair_cov = cov[index];
              values.emplace(pair_key, pair_cov);
            }

            tokens.clear();
            positions.clear();
            cov.clear();
        }
    }

    void load_samples(const string &file, vector<string>& samples) {
        ifstream input_file(file);
        string line;
        while (getline(input_file, line)) {
            samples.emplace_back(line);
        }
    }

    void flush_cache() {
        redisReply* reply = nullptr;
        reply = (redisReply*)redisCommand(redis_cache, "FLUSHDB");
        assert(reply != nullptr);
        freeReplyObject(reply);
    }
};

string LDServerTest::hostname = "127.0.0.1";
int LDServerTest::port = 6379;
boost::process::child LDServerTest::redis_server = boost::process::child();
redisContext* LDServerTest::redis_cache = nullptr;

TEST_F(LDServerTest, Morton_code) {
    ASSERT_EQ(0, to_morton_code(0, 0));
    ASSERT_EQ(42, to_morton_code(0, 7));
    ASSERT_EQ(21, to_morton_code(7, 0));
    ASSERT_EQ(63, to_morton_code(7, 7));

    uint64_t x = 0, y = 0;
    from_morton_code(63, x, y);
    ASSERT_EQ(7, x);
    ASSERT_EQ(7, y);

    from_morton_code(21, x, y);
    ASSERT_EQ(7, x);
    ASSERT_EQ(0, y);

    from_morton_code(42, x, y);
    ASSERT_EQ(0, x);
    ASSERT_EQ(7, y);

    from_morton_code(0, x, y);
    ASSERT_EQ(0, x);
    ASSERT_EQ(0, y);

    ASSERT_THROW(compute_bigmin(58, 102, 27), logic_error);
    ASSERT_EQ(compute_bigmin(58, 27, 102), 74);
    ASSERT_EQ(compute_bigmin(19, 12, 45), 36);

    ASSERT_THROW(compute_litmax(58, 102, 27), logic_error);
    ASSERT_EQ(compute_litmax(58, 27, 102), 55);
    ASSERT_EQ(compute_litmax(19, 12, 45), 15);
}

TEST_F(LDServerTest, SAV_one_page) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);

    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, correlation::LD_RSQUARE, result));
    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
    }

    result.erase();
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, correlation::LD_R, result));
    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , pow(entry.value, 2.0), 0.00000000001);
    }

    result.erase();
    //TODO: compare against pre-computed "golden standard" results. currently we jsut check if number of pairs matches the expected
    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, correlation::COV, result));
    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
//        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
    }
}

TEST_F(LDServerTest, simple_cov_test) {
    // Load gold standard covariance values from RAREMETAL
    map<string, double> gold_standard;
    this->load_raremetal_covariance("chr21.test.RAND_QT.singlevar.cov.txt", gold_standard);

    // Load score statistics (need sigma2 to multiply out)
    auto scores = load_raremetal_scores("chr21.test.RAND_QT.singlevar.score.txt");

    // LD Server start
    LDServer server(100);
    LDQueryResult result(1000);

    server.set_file("chr21.test.vcf.gz");
    ASSERT_TRUE(server.compute_region_ld("21", 9411239, 9411793, correlation::COV, result, "ALL", true));
    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        double value_gold = gold_standard.find(key)->second * scores->get_sigma2();
        double value_ldserver = entry.value;
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_NEAR(value_gold, value_ldserver, 0.0001);
    }
}

TEST_F(LDServerTest, simple_cov_test_missing_data) {
  // Load gold standard covariance values from RAREMETAL
  map<string, double> gold_standard;
  this->load_raremetal_covariance("chr21.test.missing_genotypes_and_phenotypes.RAND_QT.singlevar.cov.txt", gold_standard);

  // Load score statistics (need sigma2 to multiply out)
  auto scores = load_raremetal_scores("chr21.test.missing_genotypes_and_phenotypes.RAND_QT.singlevar.score.txt");

  // LD Server start
  LDServer server(100);
  LDQueryResult result(1000);

  server.set_file("chr21.test.missing_values.vcf.gz");
  ASSERT_TRUE(server.compute_region_ld("21", 9411239, 9411793, correlation::COV, result, "ALL", true));
  ASSERT_EQ(result.limit, 1000);
  ASSERT_EQ(result.get_last(), "");
  ASSERT_TRUE(result.is_last());
  for (auto&& entry : result.data) {
    string key(to_string(entry.position1) + "_" + to_string(entry.position2));
    double value_gold = gold_standard.find(key)->second * scores->get_sigma2();
    double value_ldserver = entry.value;
    ASSERT_NE(entry.variant1, "");
    ASSERT_NE(entry.variant2, "");
    ASSERT_NEAR(value_gold, value_ldserver, 0.005);
  }
}

TEST_F(LDServerTest, score_server) {
    // Load score statistics
    auto goldstd = load_raremetal_scores("chr21.test.RAND_QT.singlevar.score.txt");

    // Load allele frequencies
    auto goldfrq = load_variant_freqs("chr21.test.frq");

    // LD server
    LDServer ld_server(100);
    LDQueryResult ld_result(INITIAL_RESULT_SIZE);
    ld_result.limit = MAX_UINT32;
    ld_server.set_file("chr21.test.vcf.gz");

    // Score server
    ScoreServer score_server(100);
    ScoreStatQueryResult score_result(INITIAL_RESULT_SIZE);
    score_result.limit = MAX_UINT32;
    score_server.set_genotypes_file("chr21.test.vcf.gz", 1);

    ColumnTypeMap ctmap;
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    string phenotype_file = "chr21.test.tab";
    score_server.load_phenotypes_file(phenotype_file, ctmap, 2504, "\t", "iid", 1);
    score_server.set_phenotype("rand_qt");

    auto segments = make_shared_segment_vector();

    ASSERT_TRUE(ld_server.compute_region_ld("21", 9411239, 9411793, correlation::COV, ld_result, "ALL", true, segments));
    ASSERT_TRUE(score_server.compute_scores("21", 9411239, 9411793, score_result, "ALL", segments));
    ASSERT_NEAR(goldstd->get_sigma2(), score_result.sigma2, 0.001);
    for (auto&& score_res : score_result.data) {
        auto gold = goldstd->get_record(score_res.variant);
        ASSERT_NEAR(gold->u_stat, score_res.score_stat, 0.001);
        ASSERT_NEAR(gold->pvalue, score_res.pvalue, 0.001);

        auto gf = goldfrq->at(score_res.variant);
        ASSERT_NEAR(gf, score_res.alt_freq, 0.001);
    }
}

TEST_F(LDServerTest, score_server_missing_pheno) {
    // Load score statistics
    auto goldstd = load_raremetal_scores("chr21.test.missing_pheno.RAND_QT.singlevar.score.txt");

    // Load allele frequencies
    auto goldfrq = load_variant_freqs("chr21.test.frq");

    // LD server
    LDServer ld_server(100);
    LDQueryResult ld_result(INITIAL_RESULT_SIZE);
    ld_result.limit = MAX_UINT32;
    ld_server.set_file("chr21.test.vcf.gz");

    // Score server
    ScoreServer score_server(100);
    ScoreStatQueryResult score_result(INITIAL_RESULT_SIZE);
    score_result.limit = MAX_UINT32;
    score_server.set_genotypes_file("chr21.test.vcf.gz", 1);

    ColumnTypeMap ctmap;
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    string phenotype_file = "chr21.test.missing_values.tab";
    score_server.load_phenotypes_file(phenotype_file, ctmap, 2504, "\t", "iid", 1);
    score_server.set_phenotype("rand_qt");

    coordinate_samples(score_server, ld_server, "chr21.test.vcf.gz", "rand_qt", "ALL");

    auto segments = make_shared_segment_vector();

    ASSERT_TRUE(ld_server.compute_region_ld("21", 9411239, 9411793, correlation::COV, ld_result, "ALL", true, segments));
    ASSERT_TRUE(score_server.compute_scores("21", 9411239, 9411793, score_result, "ALL", segments));
    ASSERT_NEAR(goldstd->get_sigma2(), score_result.sigma2, 0.001);
    for (auto&& score_res : score_result.data) {
        auto gold = goldstd->get_record(score_res.variant);
        ASSERT_NEAR(gold->u_stat, score_res.score_stat, 0.001);
        ASSERT_NEAR(gold->pvalue, score_res.pvalue, 0.001);

        auto gf = goldfrq->at(score_res.variant);
        ASSERT_NEAR(gf, score_res.alt_freq, 0.001);
    }
}

TEST_F(LDServerTest, score_server_missing_genotypes_and_phenotypes) {
  // Load score statistics
  auto goldstd = load_raremetal_scores("chr21.test.missing_genotypes_and_phenotypes.RAND_QT.singlevar.score.txt");

  // Load allele frequencies
  auto goldfrq = load_variant_freqs("chr21.test.frq");

  // LD server
  LDServer ld_server(100);
  LDQueryResult ld_result(INITIAL_RESULT_SIZE);
  ld_result.limit = MAX_UINT32;
  ld_server.set_file("chr21.test.missing_values.vcf.gz");

  // Score server
  ScoreServer score_server(100);
  ScoreStatQueryResult score_result(INITIAL_RESULT_SIZE);
  score_result.limit = MAX_UINT32;
  score_server.set_genotypes_file("chr21.test.missing_values.vcf.gz", 1);

  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  string phenotype_file = "chr21.test.missing_values.tab";
  score_server.load_phenotypes_file(phenotype_file, ctmap, 2504, "\t", "iid", 1);
  score_server.set_phenotype("rand_qt");

  coordinate_samples(score_server, ld_server, "chr21.test.missing_values.vcf.gz", "rand_qt", "ALL");

  auto segments = make_shared_segment_vector();

  ASSERT_TRUE(ld_server.compute_region_ld("21", 9411239, 9411793, correlation::COV, ld_result, "ALL", true, segments));
  ASSERT_TRUE(score_server.compute_scores("21", 9411239, 9411793, score_result, "ALL", segments));
  ASSERT_NEAR(goldstd->get_sigma2(), score_result.sigma2, 0.001);
  for (auto&& score_res : score_result.data) {
    auto gold = goldstd->get_record(score_res.variant);
    ASSERT_NEAR(gold->u_stat, score_res.score_stat, 0.002);
    ASSERT_NEAR(gold->pvalue, score_res.pvalue, 0.002);

    auto gf = goldfrq->at(score_res.variant);
    ASSERT_NEAR(gf, score_res.alt_freq, 0.001);
  }
}

TEST_F(LDServerTest, score_covariance_runner) {
    string genotype_file = "chr22.test.vcf.gz";
    string phenotype_file = "chr22.test.tab";

    ColumnTypeMap ctmap;
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    // Region to analyze
    string chrom = "22";
    auto start = 50276998ul;
    auto stop = 50357719ul;

    // Setup mask
    vector<Mask> masks;
    Mask mask("mask.epacts.chr22.gencode-exons-AF01.tab.gz", 1, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, chrom, start, stop);
    masks.emplace_back(mask);

    // Setup ScoreCovarianceRunner configuration
    auto config = make_score_covariance_config();
    config->chrom = chrom;
    config->start = start;
    config->stop = stop;
    config->segment_size = 1000;
    config->masks = masks;
    config->sample_subset = "ALL";
    config->genotype_files = {genotype_file};
    config->genotype_dataset_id = 1;
    config->phenotype_file = phenotype_file;
    config->phenotype_column_types = ctmap;
    config->phenotype_dataset_id = 1;
    config->phenotype = "rand_qt";
    config->phenotype_nrows = 2504;
    config->phenotype_sample_column = "iid";
    config->phenotype_delim = "\t";

    // Run score/covariance calculations
    ScoreCovarianceRunner runner(config);
    runner.run();
    string json = runner.getJSON();

    // Parse back out JSON
    rapidjson::Document doc;
    doc.Parse(json.c_str());

    // Tests
    ASSERT_EQ(doc["data"]["nSamples"].GetDouble(), 2504.0);
    ASSERT_NEAR(doc["data"]["sigmaSquared"].GetDouble(), 0.08188312, 0.0001);
    ASSERT_EQ(doc["data"]["groups"][0]["variants"].Size(), 162);
    ASSERT_EQ(doc["data"]["groups"][0]["covariance"].Size(), 13203);
}

/**
 * Edge case of a group that spans a long range (2 MB) and has variants only at the ends, making many segment loads
 * unnecessary. To make it worse, we shorten the segment sizes.
 */
TEST_F(LDServerTest, scorecov_long_group_short_segments) {
  string genotype_file = "chr22.test.vcf.gz";
  string phenotype_file = "chr22.test.tab";

  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  // Region to analyze
  string chrom = "22";
  auto start = 50244251ul;
  auto stop = 51244237ul;

  // Setup mask (this would be user defined client-side via API)
  uint64_t mask_id = 0;
  vector<VariantGroup> vg_vec;
  VariantGroup vg;
  vg.chrom = "22";
  vg.name = "TEST";
  vg.start = start;
  vg.stop = stop;
  vg.variants = {
    VariantMeta("22:50244251_G/A"),
    VariantMeta("22:50244265_C/A"),
    VariantMeta("22:51244237_C/T")
  };
  vg_vec.push_back(vg);
  Mask mask(mask_id, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, vg_vec);
  vector<Mask> masks;
  masks.emplace_back(mask);

  // Setup ScoreCovarianceRunner configuration
  auto config = make_score_covariance_config();
  config->chrom = chrom;
  config->start = start;
  config->stop = stop;
  config->segment_size = 10;
  config->masks = masks;
  config->sample_subset = "ALL";
  config->genotype_files = {genotype_file};
  config->genotype_dataset_id = 1;
  config->phenotype_file = phenotype_file;
  config->phenotype_column_types = ctmap;
  config->phenotype_dataset_id = 1;
  config->phenotype = "rand_qt";
  config->phenotype_nrows = 2504;
  config->phenotype_sample_column = "iid";
  config->phenotype_delim = "\t";

  // Run score/covariance calculations
  ScoreCovarianceRunner runner(config);
  runner.run();
  string json = runner.getJSON();

  // Parse back out JSON
  rapidjson::Document doc;
  doc.Parse(json.c_str());

  // Tests
  ASSERT_EQ(doc["data"]["nSamples"].GetDouble(), 2504.0);
  ASSERT_NEAR(doc["data"]["sigmaSquared"].GetDouble(), 0.08188312, 0.0001);
  ASSERT_EQ(doc["data"]["groups"][0]["variants"].Size(), 3);
  ASSERT_EQ(doc["data"]["groups"][0]["covariance"].Size(), 6);
}

TEST_F(LDServerTest, score_user_defined_masks) {
  string genotype_file = "chr22.test.vcf.gz";
  string phenotype_file = "chr22.test.tab";

  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  // Region to analyze
  string chrom = "22";
  auto start = 50276998ul;
  auto stop = 50357719ul;

  // Setup mask (this would be user defined client-side via API)
  uint64_t mask_id = 0;
  vector<VariantGroup> vg_vec;
  VariantGroup vg;
  vg.chrom = "22";
  vg.name = "PIM3";
  vg.start = 50354416;
  vg.stop = 50357667;
  vg.variants = {VariantMeta("22:50354416_G/C"), VariantMeta("22:50355407_C/T"), VariantMeta("22:50356368_C/T"),
                 VariantMeta("22:50356386_C/T"), VariantMeta("22:50356473_C/T"), VariantMeta("22:50356497_G/A"),
                 VariantMeta("22:50356731_C/T"), VariantMeta("22:50356811_G/T"), VariantMeta("22:50356864_G/A"),
                 VariantMeta("22:50356875_C/T"), VariantMeta("22:50356887_C/T"), VariantMeta("22:50356961_C/T"),
                 VariantMeta("22:50356965_C/T"), VariantMeta("22:50356994_G/A"), VariantMeta("22:50357305_C/T"),
                 VariantMeta("22:50357350_G/A"), VariantMeta("22:50357577_G/A"), VariantMeta("22:50357657_A/G"),
                 VariantMeta("22:50357667_C/G")};
  vg_vec.push_back(vg);
  Mask mask(mask_id, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, vg_vec);
  vector<Mask> masks;
  masks.emplace_back(mask);

  // Setup ScoreCovarianceRunner configuration
  auto config = make_score_covariance_config();
  config->chrom = chrom;
  config->start = start;
  config->stop = stop;
  config->segment_size = 1000;
  config->masks = masks;
  config->sample_subset = "ALL";
  config->genotype_files = {genotype_file};
  config->genotype_dataset_id = 1;
  config->phenotype_file = phenotype_file;
  config->phenotype_column_types = ctmap;
  config->phenotype_dataset_id = 1;
  config->phenotype = "rand_qt";
  config->phenotype_nrows = 2504;
  config->phenotype_sample_column = "iid";
  config->phenotype_delim = "\t";

  // Run score/covariance calculations
  ScoreCovarianceRunner runner(config);
  runner.run();
  string json = runner.getJSON();

  // Parse back out JSON
  rapidjson::Document doc;
  doc.Parse(json.c_str());

  // Tests
  ASSERT_EQ(doc["data"]["groups"][0]["variants"].Size(), 18);
  ASSERT_EQ(doc["data"]["groups"][0]["covariance"].Size(), 171);
  ASSERT_EQ(doc["data"]["variants"][0]["variant"], "22:50354416_G/C");
  ASSERT_NEAR(doc["data"]["variants"][0]["score"].GetDouble(), -45.31790009565313, 0.0001);
  ASSERT_NEAR(doc["data"]["groups"][0]["covariance"][0].GetDouble(), 0.39843530436, 0.0001);
}

TEST_F(LDServerTest, score_no_testable_variants) {
  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  string chrom = "22";
  auto start = 50276998ul;
  auto stop = 50357719ul;

  auto config = make_score_covariance_config();
  config->chrom = chrom;
  config->start = start;
  config->stop = stop;
  config->segment_size = 1000;
  config->sample_subset = "ALL";
  config->genotype_files = {"test_no_testable_variants.vcf.gz"};
  config->genotype_dataset_id = 1;
  config->phenotype_file = "test_no_testable_variants.tab";
  config->phenotype_column_types = ctmap;
  config->phenotype_dataset_id = 1;
  config->phenotype = "rand_qt";
  config->phenotype_nrows = 2504;
  config->phenotype_sample_column = "iid";
  config->phenotype_delim = "\t";
  config->pprint();

  // Load mask
  Mask mask("test_no_testable_variants.mask.tab.gz", 1, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, chrom, start, stop);
  vector<Mask> masks;
  masks.emplace_back(mask);
  config->masks = masks;

  // Execute runner
  ScoreCovarianceRunner runner(config);
  runner.run();
  string json = runner.getJSON();

  // Parse back out JSON
  rapidjson::Document doc;
  doc.Parse(json.c_str());

  // Tests
  ASSERT_EQ(doc["data"]["groups"].Size(), 0);
  ASSERT_EQ(doc["data"]["variants"].Size(), 0);
  ASSERT_EQ(doc["data"]["nSamples"].GetDouble(), 2503.0);
}

TEST_F(LDServerTest, score_covariance_runner_monomorphic) {
    ColumnTypeMap ctmap;
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    string chrom = "22";
    auto start = 50276998ul;
    auto stop = 50357719ul;

    auto config = make_score_covariance_config();
    config->chrom = chrom;
    config->start = start;
    config->stop = stop;
    config->segment_size = 1000;
    config->sample_subset = "ALL";
    config->genotype_files = {"chr22.monomorphic_test.vcf.gz"};
    config->genotype_dataset_id = 1;
    config->phenotype_file = "chr22.test.missing_values.tab";
    config->phenotype_column_types = ctmap;
    config->phenotype_dataset_id = 1;
    config->phenotype = "rand_qt";
    config->phenotype_nrows = 2504;
    config->phenotype_sample_column = "iid";
    config->phenotype_delim = "\t";
    config->pprint();

    // Load mask
    Mask mask("mask.epacts.chr22.gencode-exons-AF01.tab.gz", 1, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, chrom, start, stop);
    vector<Mask> masks;
    masks.emplace_back(mask);
    config->masks = masks;

    // Execute runner
    ScoreCovarianceRunner runner(config);
    runner.run();
    string json = runner.getJSON();

    // Parse back out JSON
    rapidjson::Document doc;
    doc.Parse(json.c_str());

    // Tests
    ASSERT_EQ(doc["data"]["nSamples"].GetDouble(), 2495.0);
    ASSERT_NEAR(doc["data"]["sigmaSquared"].GetDouble(), 0.081971, 0.0001);
    ASSERT_EQ(doc["data"]["variants"].Size(), 685);
    ASSERT_EQ(doc["data"]["groups"][0]["variants"].Size(), 159);
    ASSERT_EQ(doc["data"]["groups"][0]["covariance"].Size(), 12720);

    auto& variants = doc["data"]["variants"];
    set<string> all_variants;

    for (auto&& vblock : variants.GetArray()) {
      string variant = vblock.GetObject()["variant"].GetString();
      all_variants.emplace(variant);
    }

    auto& groups = doc["data"]["groups"];
    for (auto&& group : groups.GetArray()) {
      uint64_t n_variants = group.GetObject()["variants"].Size();
      uint64_t n_covar = group.GetObject()["covariance"].Size();
      ASSERT_EQ(n_covar, n_variants * (n_variants + 1) / 2);

      auto& group_variants = group.GetObject()["variants"];
      for (auto&& vobj : group_variants.GetArray()) {
        string variant = vobj.GetString();
        ASSERT_TRUE(all_variants.find(variant) != all_variants.end());
      }
    }
}

TEST_F(LDServerTest, BCF_one_page) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);

    server.set_file("chr22.test.bcf");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, correlation::LD_RSQUARE, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
    }
}


TEST_F(LDServerTest, VCF_one_page) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);

    server.set_file("chr22.test.vcf.gz");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, correlation::LD_RSQUARE, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
    }
}

TEST_F(LDServerTest, SAV_chrX_one_page) {
    map<string, double> goldstandard;
    vector<string> samples;
    this->load_region_goldstandard("region_ld_X_60100_60150.hap.ld", goldstandard);
    this->load_samples("EUR.samples.txt", samples);

    LDServer server(100);
    LDQueryResult result(1000);

    server.set_file("chrX.test.sav");
    ASSERT_TRUE(server.compute_region_ld("X", 60100, 60150, correlation::LD_RSQUARE, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.000000001);
    }
}

TEST_F(LDServerTest, region_with_paging) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(4);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_region_ld("22", 51241101, 51241385, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 4);
        ASSERT_LE(result.data.size(), 4);
        for (auto &&entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, variant_with_paging_1) {
    map<string, double> goldstandard;
    this->load_variant_goldstandard("variant_ld_22_51241101_vs_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(2);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:51241101_A/T", "22", 51241101, 51241385, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
        for (auto&& entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, variant_with_paging_2) {
    map<string, double> goldstandard;
    this->load_variant_goldstandard("variant_ld_22_51241386_vs_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(2);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:51241386_C/G", "22", 51241101, 51241385, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
        for (auto&& entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.0000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, variant_with_paging_3) {
    map<string, double> goldstandard;
    this->load_variant_goldstandard("variant_ld_22_51241309_vs_51241101_51244237.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(2);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:51241309_C/T", "22", 51241101, 51244237, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
        for (auto&& entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.0000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, variant_with_paging_no_variant_1) {
    LDServer server(100);
    LDQueryResult result(2);
    int result_total_size = 0;
    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:51241101_A/C", "22", 51241101, 51241385, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, 0);
}

TEST_F(LDServerTest, variant_with_paging_no_variant_2) {
    LDServer server(100);
    LDQueryResult result(2);
    int result_total_size = 0;
    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:51241386_C/T", "22", 51241101, 51241385, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, 0);
}

TEST_F(LDServerTest, variant_with_paging_no_variant_3) {
    LDServer server(100);
    LDQueryResult result(2);
    int result_total_size = 0;
    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:51241309_C/A", "22", 51241101, 51244237, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
    }
    ASSERT_EQ(result_total_size, 0);
}


TEST_F(LDServerTest, AFR_region_with_paging) {
    map<string, double> goldstandard;
    vector<string> samples;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.AFR.hap.ld", goldstandard);
    this->load_samples("AFR.samples.txt", samples);

    LDServer server(100);
    LDQueryResult result(2);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    server.set_samples("AFR", samples);
    while (server.compute_region_ld("22", 51241101, 51241385, correlation::LD_RSQUARE, result, "AFR")) {
        ASSERT_LE(result.limit, 4);
        ASSERT_LE(result.data.size(), 4);
        for (auto &&entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.0000000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, large_region_with_paging) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_50544251_50549251.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_region_ld("22", 50544251, 50549251, correlation::LD_RSQUARE, result, "ALL")) {
        ASSERT_LE(result.limit, 1000);
        ASSERT_LE(result.data.size(), 1000);
        for (auto &&entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, large_variant_with_paging) {
    map<string, double> goldstandard;
    this->load_variant_goldstandard("variant_ld_22_50546666_vs_50544251_50549251.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:50546666_C/T", "22", 50544251, 50549251, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 1000);
        ASSERT_LE(result.data.size(), 1000);
        for (auto&& entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, hiredis) {
    redisContext *context = nullptr;
    struct timeval timeout = { 1, 500000 }; // 1.5 seconds
    context = redisConnectWithTimeout(hostname.c_str(), port, timeout);
    ASSERT_NE(context, nullptr);
    ASSERT_EQ(context->err, 0);
    redisFree(context);
}

TEST_F(LDServerTest, segment_key) {
    uint32_t unique_key = 0u;
    uint64_t start_bp = 0u, stop_bp = 0u;

    string key = LDServer::make_segment_cache_key(1, "", "", 0, 1);
    ASSERT_EQ(20, key.size());
    memcpy(&unique_key, key.c_str(), 4);
    memcpy(&start_bp, key.c_str() + 4, 8);
    memcpy(&stop_bp, key.c_str() + 12, 8);
    ASSERT_EQ(unique_key, 1);
    ASSERT_EQ(start_bp, 0);
    ASSERT_EQ(stop_bp, 1);

    key = LDServer::make_segment_cache_key(2, "ALL", "chr22", 10, 20);
    ASSERT_EQ(28, key.size());
    memcpy(&unique_key, key.c_str(), 4);
    ASSERT_EQ(unique_key, 2);
    string samples_name(key.c_str() + 4, 3);
    ASSERT_EQ(samples_name, "ALL");
    string chromosome(key.c_str() + 7, 5);
    ASSERT_EQ(chromosome, "chr22");
    memcpy(&start_bp, key.c_str() + 12 , 8);
    memcpy(&stop_bp, key.c_str() + 20, 8);
    ASSERT_EQ(start_bp, 10);
    ASSERT_EQ(stop_bp, 20);
}

TEST_F(LDServerTest, segment_cache) {
    RawSAV raw("chr22.test.sav");
    raw.open("22", raw.get_samples(), false);
    string key = LDServer::make_segment_cache_key(1, "ALL", "22", 51241101, 51241385);
    Segment segment1("22", 51241101, 51241385, genotypes_store::CSC_ALL_ONES);
    ASSERT_FALSE(segment1.has_names());
    ASSERT_FALSE(segment1.has_genotypes());
    raw.load(segment1);
    ASSERT_TRUE(segment1.has_names());
    ASSERT_TRUE(segment1.has_genotypes());

    segment1.save(redis_cache, key);

    Segment segment2("22", 51241101, 51241385, genotypes_store::CSC_ALL_ONES);
    ASSERT_FALSE(segment2.has_names());
    ASSERT_FALSE(segment2.has_genotypes());
    ASSERT_EQ(segment2.get_n_variants(), 0);
    ASSERT_EQ(segment2.get_n_haplotypes(), 0);
    segment2.load(redis_cache, key);
    ASSERT_TRUE(segment2.has_names());
    ASSERT_FALSE(segment2.has_genotypes());
    ASSERT_EQ(segment2.get_n_variants(), segment1.get_n_variants());
    for (int i = 0; i < segment2.get_n_variants(); ++i) {
        ASSERT_EQ(segment2.get_name(i), segment1.get_name(i));
        ASSERT_EQ(segment2.get_position(i), segment1.get_position(i));
    }
}

TEST_F(LDServerTest, cell_key) {
    uint32_t unique_key = 0;
    correlation correlation_type = correlation::LD_R;
    uint64_t morton_code = 0;

    string key = LDServer::make_cell_cache_key(1, "", correlation::LD_RSQUARE, "", 3);

    ASSERT_EQ(13, key.size());
    memcpy(&unique_key, key.c_str(), 4);
    ASSERT_EQ(unique_key, 1);
    memcpy(&correlation_type, key.c_str() + 4, 1);
    ASSERT_EQ(correlation_type, correlation::LD_RSQUARE);
    memcpy(&morton_code, key.c_str() + 5, 8);
    ASSERT_EQ(morton_code, 3);

    key = LDServer::make_cell_cache_key(2, "ALL", correlation::COV, "chr22", 300);
    ASSERT_EQ(21, key.size());
    memcpy(&unique_key, key.c_str(), 4);
    ASSERT_EQ(unique_key, 2);
    string samples_name(key.c_str() + 4, 3);
    ASSERT_EQ(samples_name, "ALL");
    string chromosome(key.c_str() + 7, 5);
    ASSERT_EQ(chromosome, "chr22");
    memcpy(&correlation_type, key.c_str() + 12, 1);
    ASSERT_EQ(correlation_type, correlation::COV);
    memcpy(&morton_code, key.c_str() + 13 , 8);
    ASSERT_EQ(morton_code, 300);
}

TEST_F(LDServerTest, cell_cache) {
    LDQueryResult result1(1000);
    LDQueryResult result2(1000);

    RawSAV raw("chr22.test.sav");
    raw.open("22", raw.get_samples(), false);
    string cell_key = LDServer::make_cell_cache_key(1, "ALL", correlation::LD_RSQUARE, "22", to_morton_code(512411, 512411));
    string segment_key = LDServer::make_segment_cache_key(1, "ALL", "22", 51241100, 51241199);

    CellFactory factory;

    shared_ptr<Cell> cell1 = factory.create(correlation::LD_RSQUARE, 512411, 512411);
    cell1->segment_i = make_shared<Segment>("22", 51241100, 51241199, genotypes_store::CSC_ALL_ONES);
    raw.load(*(cell1->segment_i));
    cell1->segment_i->save(redis_cache, segment_key);
    cell1->compute();
    cell1->save(redis_cache, cell_key);
    cell1->extract(51241100, 51241199, result1);

    shared_ptr<Cell> cell2 = factory.create(correlation::LD_RSQUARE, 512411, 512411);
    cell2->segment_i = make_shared<Segment>("22", 51241100, 51241199, genotypes_store::CSC_ALL_ONES);
    cell2->segment_i->load(redis_cache, segment_key);
    cell2->load(redis_cache, cell_key);
    cell2->extract(51241100, 51241199, result2);

    ASSERT_EQ(cell1->segment_i->get_n_variants(), cell2->segment_i->get_n_variants());
    for (unsigned int i = 0; i < cell1->segment_i->get_n_variants(); ++i) {
        ASSERT_EQ(cell1->segment_i->get_name(i), cell2->segment_i->get_name(i));
        ASSERT_EQ(cell1->segment_i->get_position(i), cell2->segment_i->get_position(i));
    }

    ASSERT_EQ(result1.data.size(), result2.data.size());
    for (unsigned int i = 0; i < result1.data.size(); ++i) {
        ASSERT_EQ(result1.data[i].variant1, result2.data[i].variant1);
        ASSERT_EQ(result1.data[i].variant2, result2.data[i].variant2);
        ASSERT_EQ(result1.data[i].value, result2.data[i].value);
    }
}

TEST_F(LDServerTest, cache_enabled) {
    redisReply* reply = nullptr;
    LDServer server;
    LDQueryResult result(1000);

    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241199, correlation::LD_RSQUARE, result));

    reply = (redisReply*)redisCommand(redis_cache, "DBSIZE");
    ASSERT_NE(reply, nullptr);
    ASSERT_EQ(reply->integer, 0);
    freeReplyObject(reply);

    result.erase();
    server.enable_cache(1, hostname, port);
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241199, correlation::LD_RSQUARE, result));
    reply = (redisReply*)redisCommand(redis_cache, "DBSIZE");
    ASSERT_NE(reply, nullptr);
    ASSERT_EQ(reply->integer, 2);
    freeReplyObject(reply);
}

TEST_F(LDServerTest, cache_with_different_correlations) {
    LDServer server;
    LDQueryResult result1(1000);
    LDQueryResult result2(1000);

    flush_cache();

    server.enable_cache(1, hostname, port);
    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 50248762, 50249221, correlation::LD_RSQUARE, result1));
    ASSERT_TRUE(server.compute_region_ld("22", 50248762, 50249221, correlation::LD_R, result2));
    ASSERT_TRUE(result1.is_last());
    ASSERT_TRUE(result2.is_last());
    ASSERT_EQ(result1.data.size(), result2.data.size());
    for (unsigned int i = 0u; i < result1.data.size(); ++i) {
        ASSERT_NE(result1.data[i].value, result2.data[i].value);
    }

    flush_cache();
}

TEST_F(LDServerTest, large_region_with_paging_and_caching) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_50544251_50549251.hap.ld", goldstandard);

    LDServer server(1000);
    server.enable_cache(1, hostname.c_str(), port);
    LDQueryResult result(1000);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_region_ld("22", 50544251, 50549251, correlation::LD_RSQUARE, result, "ALL")) {
        ASSERT_LE(result.limit, 1000);
        ASSERT_LE(result.data.size(), 1000);
        for (auto &&entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, SAV_one_page_no_segment_intersect_1) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);

    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241280, correlation::LD_RSQUARE, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), 1);
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
    }
}

TEST_F(LDServerTest, SAV_one_page_no_segment_intersect_2) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);

    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 51241103, 51241385, correlation::LD_RSQUARE, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_GT(result.data.size(), 1);
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
    }
}

TEST_F(LDServerTest, known_memory_leak_test) {
    LDServer server(1000);
    server.enable_cache(1, hostname.c_str(), port);
    LDQueryResult result(100000);
    server.set_file("../../data/ALL.chr22.sav");
    while (server.compute_variant_ld("22:45186606_C/T", "22", 45186270, 45286270, correlation::LD_RSQUARE, result, LDServer::ALL_SAMPLES_KEY));
}

TEST_F(LDServerTest, empty_data_test) {
    LDServer server(100);
    LDQueryResult result(1000);
    server.set_file("chr22.test.sav");
    ASSERT_FALSE(server.compute_region_ld("21", 51241103, 51241385, correlation::LD_RSQUARE, result));
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), 0);
}

TEST_F(LDServerTest, rsquare_approx_test) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server(100);
    LDQueryResult result(1000);
    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, correlation::LD_RSQUARE_APPROX, result));
    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_TRUE(result.is_last());
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.value, 0.00000000001);
    }
}

TEST_F(LDServerTest, DISABLED_read_spead) {
    vector<string> names;
    vector<uint64_t> positions;

    vector<arma::uword> sp_mat_rowind;
    vector<arma::uword> sp_mat_colind;

    uint64_t n_haplotypes;

    auto start = std::chrono::system_clock::now();
    savvy::indexed_reader f("../../data/ALL.chr22.sav", {"22", 25637000, 26137999}, savvy::fmt::gt);
    std::chrono::duration<double> elapsed =  std::chrono::system_clock::now() - start;
    std::cout << "reader setup: " << elapsed.count() << " s\n";


    std::stringstream ss;
    savvy::site_info anno;
    savvy::armadillo::sparse_vector<float> alleles;
    while (f.read(anno, alleles)) {
        n_haplotypes = alleles.n_rows;
        if (alleles.n_nonzero > 0) {
            ss.str("");
            ss << anno.chromosome() << ":" << anno.position() << "_" << anno.ref() << "/" << anno.alt();
            names.emplace_back(ss.str());
            positions.push_back(anno.position());
            sp_mat_colind.push_back(sp_mat_rowind.size());
            for (auto it = alleles.begin(); it != alleles.end(); ++it) {
                sp_mat_rowind.push_back(it.row());
            }
        }
    }
    sp_mat_colind.push_back(sp_mat_rowind.size());

    cout << names.size() << endl;
    cout << positions.size() << endl;
}

TEST_F(LDServerTest, pheno_read_tab) {
    Phenotypes pheno;
    ColumnTypeMap ctmap;
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    pheno.load_file("chr21.test.tab", ctmap, 2504, "\t", "iid");

    auto &rand_qt = *pheno.as_vec("rand_qt");
    ASSERT_NEAR(rand_qt[0], 0.283211535632, 0.000001);

    auto &rand_binary = *pheno.as_vec("rand_binary");
    ASSERT_TRUE(isnan(rand_binary[0]));
    ASSERT_EQ(rand_binary[1], 1);
    ASSERT_EQ(rand_binary[2], 0);

    auto &sex = *pheno.as_vec("sex");
    ASSERT_EQ(sex[0], 0);
    ASSERT_EQ(sex[1], 1);

    auto phenos_loaded = *pheno.get_phenotypes();
    ASSERT_EQ(phenos_loaded[0], "iid");
}

TEST_F(LDServerTest, pheno_for_analysis_columns) {
  Phenotypes pheno;
  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::FLOAT); // if this column is parsed, it will fail b/c set to float when it is "male"/"female"
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  auto analysis_cols = make_shared<vector<string>>(initializer_list<string>{"rand_qt", "rand_binary"});
  pheno.load_file("chr22.test.tab", ctmap, 2504, "\t", "iid", analysis_cols);
  auto phenotypes = pheno.get_phenotypes();
  ASSERT_EQ(phenotypes->size(), 3);
}

TEST_F(LDServerTest, pheno_bad_float) {
  using ::testing::HasSubstr;

  Phenotypes pheno;
  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  try {
    pheno.load_file("chr22.test.bad_float.tab", ctmap, 2504, "\t", "iid");
    FAIL() << "Expected std::runtime_error";
  }
  catch (PhenotypeParseException& e) {
    EXPECT_THAT(e.what(), HasSubstr("Error reading line 3, column 3 (rand_qt) in phenotype file"));
    EXPECT_THAT(e.what(), HasSubstr("invalid value: something"));
  }
  catch (...) {
    FAIL() << "Expected std::runtime_error";
  }
}

TEST_F(LDServerTest, pheno_read_ped) {
    Phenotypes pheno;
    ColumnTypeMap ctmap;

    ctmap.add("fid", ColumnType::TEXT);
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("patid", ColumnType::TEXT);
    ctmap.add("matid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    pheno.load_file("chr21.test.ped", ctmap, 2504, "\t", "iid");

    auto &rand_qt = *pheno.as_vec("rand_qt");
    ASSERT_NEAR(rand_qt[0], 0.283211535632, 0.000001);

    auto &rand_binary = *pheno.as_vec("rand_binary");
    ASSERT_TRUE(isnan(rand_binary[0]));
    ASSERT_EQ(rand_binary[1], 1);
    ASSERT_EQ(rand_binary[2], 0);

    auto &sex = *pheno.as_vec("sex");
    ASSERT_EQ(sex[0], 1);
    ASSERT_EQ(sex[1], 0);

    auto phenos_loaded = *pheno.get_phenotypes();
    ASSERT_EQ(phenos_loaded[0], "fid");
}

TEST_F(LDServerTest, pheno_reorder) {
    Phenotypes pheno;
    ColumnTypeMap ctmap;
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    pheno.load_file("chr21.test.tab", ctmap, 2504, "\t", "iid");

    // Reorder samples
    vector<string> new_samples = {"HG00100","HG00103","HG00096","BAD_SAMPLE"};
    pheno.reorder(new_samples);

    // Check values
    auto &rand_qt = *pheno.as_vec("rand_qt");
    ASSERT_NEAR(rand_qt[0], 0.2846738, 0.000001);
    ASSERT_NEAR(rand_qt[1], 0.7341857, 0.000001);
    ASSERT_NEAR(rand_qt[2], 0.2832115, 0.000001);
    ASSERT_TRUE(isnan(rand_qt[3]));
}

/**
 * Currently this is just a quick regression test, not verified for correctness.¡
 */
TEST_F(LDServerTest, pheno_compute_score) {
    Phenotypes pheno;
    ColumnTypeMap ctmap;
    ctmap.add("iid", ColumnType::TEXT);
    ctmap.add("sex", ColumnType::CATEGORICAL);
    ctmap.add("rand_binary", ColumnType::CATEGORICAL);
    ctmap.add("rand_qt", ColumnType::FLOAT);

    pheno.load_file("chr21.test.tab", ctmap, 2504, "\t", "iid");

    vector<double> test_geno = {0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,1,0,0,1,0,0,0,0,2,0,1,0,
                               1,0,0,0,1,1,1,0,1,0,0,1,1,1,0,1,0,1,0,0,0,0,1,0,0,1,0,0,2,1,2,1,0,0,0,
                               0,1,1,1,1,0,0,1,1,0,0,1,0,0,1,0,0,0,0,1,2,1,0,0,0,0,1,0,1,0,1,0,0,0,0,
                               1,0,1,0,1,1,1,0,1,0,1,0,0,1,0,0,2,0,0,1,0,1,0,0,1,0,1,0,1,1,0,1,1,0,1,
                               1,0,0,0,0,2,0,0,1,1,1,0,1,0,1,0,0,1,1,2,0,1,1,1,0,1,1,0,1,0,0,0,0,0,1,
                               0,0,1,0,0,0,0,0,1,1,0,0,2,1,0,1,1,0,1,1,0,0,0,0,0,1,1,1,0,1,0,2,1,1,1,
                               1,1,0,0,1,1,0,1,1,1,2,1,1,0,1,0,0,2,0,0,1,0,0,1,1,0,0,0,1,0,0,1,0,1,1,
                               2,0,2,0,0,1,0,0,0,1,1,1,1,0,1,1,1,2,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,
                               0,1,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,2,0,1,0,0,0,0,1,0,
                               0,1,1,1,0,2,1,1,0,1,0,0,0,0,0,0,1,0,1,2,0,2,1,0,0,0,0,0,0,1,1,1,0,1,0,
                               2,0,0,2,0,0,0,0,1,2,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,2,1,0,1,
                               0,0,0,1,1,0,0,0,0,1,0,1,0,0,1,0,2,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,0,
                               1,0,1,0,0,0,1,0,0,1,1,0,2,0,0,0,1,0,1,1,0,1,0,0,1,1,1,1,0,0,0,0,0,1,0,
                               0,0,1,1,0,1,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,2,0,0,1,2,1,
                               0,0,0,1,1,0,2,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,0,0,
                               1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,2,2,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,
                               1,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,1,1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,2,
                               0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,2,1,1,1,0,0,0,0,0,0,0,
                               0,1,1,1,2,0,0,1,1,2,0,0,0,1,0,1,1,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,2,
                               1,0,2,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,1,2,0,1,0,0,1,1,1,0,0,1,1,1,1,
                               2,1,0,0,0,1,1,0,0,0,2,1,1,0,0,0,1,0,0,1,1,1,1,0,1,0,0,1,1,1,0,1,0,0,0,
                               0,0,1,1,1,0,1,1,0,1,0,0,1,1,1,1,2,1,0,0,1,0,2,1,0,1,2,1,0,0,2,0,0,1,1,
                               1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,2,1,1,1,2,0,0,0,0,0,0,0,0,2,1,0,1,1,
                               0,0,1,0,1,0,2,0,0,0,1,0,1,0,0,1,1,0,1,1,0,0,1,2,0,0,0,0,0,0,0,1,0,0,1,
                               0,1,1,1,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,2,1,1,0,0,0,1,2,0,0,0,1,1,0,
                               1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,0,0,0,
                               1,0,0,0,0,0,0,1,0,0,0,0,2,1,1,1,1,1,1,1,0,0,2,0,1,0,1,1,1,0,0,0,1,1,0,
                               1,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,0,0,0,2,0,1,0,0,0,1,1,0,0,0,1,1,2,
                               2,2,0,0,0,2,1,1,0,0,0,1,1,0,1,1,0,1,1,1,1,0,2,2,1,1,0,1,1,2,1,1,0,1,0,
                               1,0,1,0,1,0,1,0,0,1,1,0,0,1,0,0,0,1,2,0,1,0,1,0,1,1,1,0,1,1,0,0,0,0,0,
                               1,1,0,0,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,2,0,1,0,0,1,
                               0,1,0,0,0,0,1,1,0,1,0,1,2,0,0,0,0,0,2,0,0,0,1,2,0,0,0,1,1,0,1,1,0,1,1,
                               0,0,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,0,1,1,0,0,0,1,0,0,1,0,0,1,0,0,1,1,2,
                               1,1,0,1,1,0,1,1,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,0,
                               0,1,1,0,0,1,0,2,0,0,0,1,0,0,0,0,2,1,1,1,0,0,0,2,0,1,0,0,2,0,0,1,0,0,2,
                               0,1,1,1,0,1,0,0,0,0,0,2,1,0,1,0,1,1,0,0,0,0,2,1,2,1,1,0,1,1,0,2,0,0,2,
                               0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,1,0,0,0,1,0,0,1,1,1,0,1,2,0,2,1,1,0,1,2,
                               0,0,0,1,0,0,1,1,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,0,1,1,0,0,1,0,2,1,
                               1,1,1,0,1,1,0,0,1,1,0,2,1,1,0,0,2,1,0,1,0,0,0,0,1,1,0,1,1,1,0,1,0,1,1,
                               0,0,1,1,1,1,1,2,0,2,1,0,0,0,2,1,0,1,1,1,1,0,1,0,0,0,1,0,0,0,1,0,0,0,1,
                               0,0,0,0,0,1,0,0,0,0,2,0,1,0,0,2,0,1,0,1,0,1,0,1,0,0,2,0,0,1,0,1,0,0,0,
                               1,0,2,0,0,0,1,0,1,0,2,0,0,0,2,1,0,0,0,2,1,1,0,0,1,1,0,0,2,1,0,0,0,1,1,
                               0,0,0,0,0,0,1,0,1,1,0,1,0,0,0,1,2,1,0,1,0,0,0,0,2,0,1,1,1,2,2,1,1,0,1,
                               0,1,1,1,1,1,1,2,1,0,1,1,0,0,0,1,0,1,1,2,2,0,1,1,0,0,0,0,1,0,1,1,0,0,0,
                               0,0,0,1,0,1,0,0,0,1,1,1,1,0,0,0,0,1,1,0,1,0,0,1,0,1,1,0,0,0,1,0,1,1,0,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0,0,1,1,1,1,0,0,1,1,
                               0,0,0,1,2,0,0,1,0,1,2,0,0,2,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,2,1,0,
                               1,0,0,0,0,0,1,1,1,0,0,2,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,1,0,2,0,0,0,1,
                               1,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,2,0,1,0,
                               0,1,1,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,0,1,2,0,0,1,0,0,1,0,0,1,1,1,0,1,1,
                               1,0,1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,0,0,1,0,0,1,1,1,
                               0,1,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1,0,2,1,2,1,0,0,2,0,0,0,
                               0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,2,0,2,1,
                               1,0,2,1,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0,1,0,0,0,1,1,0,2,0,2,0,0,0,0,1,
                               0,0,1,0,0,2,1,1,0,1,1,1,0,2,0,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,1,1,1,2,0,
                               0,0,2,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,1,1,0,1,1,1,
                               0,0,0,1,0,1,1,1,2,2,1,1,1,0,1,1,0,1,1,0,0,0,1,1,0,0,0,2,0,0,0,0,1,0,1,
                               0,1,0,0,1,0,2,0,0,1,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,1,0,0,0,1,
                               0,0,0,0,1,0,1,0,0,1,1,1,2,0,0,0,2,1,1,2,0,1,0,2,0,0,1,0,1,1,1,0,0,1,2,
                               0,1,0,0,0,1,1,0,0,1,1,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,2,2,1,1,
                               0,0,0,0,1,1,0,0,0,2,1,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,2,1,0,0,1,1,1,0,1,
                               0,1,0,0,0,0,1,1,1,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,2,
                               2,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,1,0,2,1,
                               1,0,1,0,0,2,2,0,0,0,0,1,0,1,0,2,1,0,1,1,1,0,1,1,0,1,0,1,0,1,0,2,1,0,0,
                               1,1,0,1,1,0,0,0,0,0,0,1,0,1,1,1,0,1,1,0,0,1,0,1,0,0,0,2,0,1,0,0,0,1,0,
                               0,2,1,1,0,1,1,0,0,0,0,0,1,0,2,0,0,0,0,2,0,0,1,0,0,0,1,1,0,0,1,0,0,0,1,
                               0,0,1,1,0,1,1,0,1,0,2,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,1,0,1,0,0,
                               0,0,1,1,1,0,1,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,1,0,2,0,0,0,1,
                               0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,0,1,0,2,0,1,2,0,1,1,1,1,1,0,0,1,1,1,0,
                               0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,1,0,1,2,1,0,1,1,2,0,0,0,
                               1,0,1,1,0,0,2,0,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,2,1,1,0,0,1,
                               0,0,0,0,0,0,1,1,0,0,0,1,0,0,2,1,1,0,0};

    arma::vec geno_vec(test_geno);
    auto result = pheno.compute_score(geno_vec, "rand_qt");
    ASSERT_NEAR(result->score_stat, 0.37527933, 0.000001);
    ASSERT_NEAR(pheno.compute_sigma2("rand_qt"), 0.08188312, 0.000001);
    ASSERT_NEAR(result->pvalue, 0.18969995, 0.000001);
}

TEST_F(LDServerTest, score_server_paging) {
  LDServer ld_server(100);
  LDQueryResult ld_result(1000);

  ScoreServer score_server(100);
  ScoreStatQueryResult score_results(3);
  auto segments = make_shared_segment_vector();

  ld_server.set_file("chr22.test.vcf.gz");
  ld_server.compute_region_ld("22", 51241101, 51241385, correlation::COV, ld_result, "ALL", true, segments);

  score_server.set_genotypes_file("chr22.test.vcf.gz", 1);

  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);
  score_server.load_phenotypes_file("chr21.test.tab", ctmap, 2504, "\t", "iid", 1);
  score_server.set_phenotype("rand_qt");
  score_server.compute_scores("22", 51241101, 51241385, score_results, "ALL", segments);

  ASSERT_EQ(score_results.data.size(),3);
  ASSERT_EQ(score_results.last_seg,1);
  ASSERT_EQ(score_results.last_i,1);
  ASSERT_EQ(score_results.page,1);
}

TEST_F(LDServerTest, score_segment_cache) {
    RawSAV raw("chr22.test.sav");
    raw.open("22", raw.get_samples(), false);

    string key = ScoreServer::make_segment_cache_key(1, 1, "rand_qt", "ALL", "22", 51241101, 51241385);

    ScoreSegment segment1("22", 51241101, 51241385, genotypes_store::CSC);
    ASSERT_FALSE(segment1.has_names());
    ASSERT_FALSE(segment1.has_genotypes());

    raw.load(segment1);
    ASSERT_TRUE(segment1.has_names());
    ASSERT_TRUE(segment1.has_genotypes());

    ScoreResult result;
    result.variant = "22:51241101";
    result.pvalue = 0.454;
    result.score_stat = 0.454;
    result.alt_freq = 0.454;

    segment1.add_score(result);

    segment1.save(redis_cache, key);

    ScoreSegment segment2("22", 51241101, 51241385, genotypes_store::CSC);
    ASSERT_FALSE(segment2.has_names());
    ASSERT_FALSE(segment2.has_genotypes());
    ASSERT_EQ(segment2.get_n_variants(), 0);
    ASSERT_EQ(segment2.get_n_haplotypes(), 0);

    segment2.load(redis_cache, key);

    ASSERT_TRUE(segment2.has_names());
    ASSERT_FALSE(segment2.has_genotypes());
    ASSERT_EQ(segment2.get_n_variants(), segment1.get_n_variants());
    for (int i = 0; i < segment2.get_n_variants(); ++i) {
        ASSERT_EQ(segment2.get_name(i), segment1.get_name(i));
        ASSERT_EQ(segment2.get_position(i), segment1.get_position(i));
    }

    ASSERT_TRUE(segment2.has_scores());
    ASSERT_TRUE(segment1 == segment2);
}

TEST_F(LDServerTest, regression_logistic_small) {
  arma::mat x = {
    { 1, -0.096804537415269, 0.394521501456981, 0.136719199203261  },
    { 1, 1.87480075412442, -0.576503014867386, -0.294691393024923  },
    { 1, 0.772754744319657, -0.558625773661303, -0.976660217651494  },
    { 1, 0.717732517365932, 0.0936844055088714, 1.78662210708278  },
    { 1, 1.8018548122351, 1.85515254706929, -1.70359905461682  },
    { 1, -0.348127414820295, -3.05083075001173, -0.0369152710070338  },
    { 1, -1.01502457917055, -0.276941129747707, -0.406115851940272  },
    { 1, 1.10078194115342, 0.515343122857823, 0.995075331192894  },
    { 1, 0.34158930036517, 0.570973367742649, 0.817951168123065  },
    { 1, 0.778222729194391, -0.717535882710764, 0.276220595509114  },
    { 1, -0.445580619793771, -1.35355138636549, 0.0400984960691716  },
    { 1, 1.05968330903572, 2.48073248812764, 0.744166070387559  },
    { 1, 1.13193063729192, 0.803345354713477, 0.67083197357916  },
    { 1, 0.64438005748593, -0.820235260483927, -1.52106595803099  },
    { 1, 0.998139980346589, -0.101007365065264, 0.742100369064892  },
    { 1, 0.0043382638725237, 1.44109213606475, -0.0368974515646374  },
    { 1, 1.2213571153, -0.547334882194856, -0.947655460648014  },
    { 1, -0.545318955493577, -0.167697594444228, 1.18861767550958  },
    { 1, -0.11332905303665, 0.477938596320964, -0.124355176118656  },
    { 1, -1.06063645191269, -2.49783021515651, -0.579642166533558  },
    { 1, -1.31296784040276, -0.597392049046657, 0.0267899494002477  },
    { 1, -1.14312422332763, 2.12767104213113, -0.223249127365326  },
    { 1, 0.832580243798403, -1.78912384982458, -1.00652278049377  },
    { 1, -1.12411740497594, -0.795110424818677, -0.947332980605162  },
    { 1, 0.705237914919055, 0.222294866312317, -0.39563859497371  },
    { 1, 0.229680544111858, 1.25628142318361, 0.741735315848646  },
    { 1, -1.0101433303265, -1.23653274314961, -0.291369194678221  },
    { 1, 0.82667455981141, 1.33750612592533, -1.46654503496563  },
    { 1, 0.528798998554926, 0.190260018116963, -1.13765940889545  },
    { 1, 1.18406942883136, 0.851258176994322, -1.74868221717553  },
  };

  arma::vec y = { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0};
  LogisticRegression model;
  model.fit(y, x);

  /*
    Call:
    glm(formula = y ~ x, family = "binomial")

    Deviance Residuals:
         Min        1Q    Median        3Q       Max
    -1.54734  -1.10953   0.06771   1.08774   1.60218

    Coefficients:
                Estimate Std. Error z value Pr(>|z|)
    (Intercept) -0.18626    0.42311  -0.440   0.6598
    x1           0.55630    0.49250   1.130   0.2587
    x2          -0.67281    0.38355  -1.754   0.0794 .
    x3          -0.07457    0.43093  -0.173   0.8626
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

    (Dispersion parameter for binomial family taken to be 1)

        Null deviance: 41.589  on 29  degrees of freedom
    Residual deviance: 37.319  on 26  degrees of freedom
    AIC: 45.319

    Number of Fisher Scoring iterations: 4
   */

  auto betas = *model.getBetas();
  auto std_err = *model.getStandardErrors();
  auto p = *model.getPredictedProb();
  auto resid = *model.getResiduals();
  auto pvals = *model.getPvalues();

  vector<double> expected_p = { 0.373846362774307, 0.780149601877233, 0.666478653815652, 0.504185638464105, 0.424329673684531, 0.842303173334388, 0.369513437699039, 0.50129076353375, 0.391408829293206, 0.670142273514121, 0.616225441901857, 0.210597307527924, 0.463290963169013, 0.697938227509708, 0.594263151190782, 0.240364041411177, 0.71750074604553, 0.385697671445548, 0.363177706176422, 0.72061319902535, 0.373617316780057, 0.0964716804290593, 0.82573608136868, 0.448680997567444, 0.521486628711504, 0.277068575402555, 0.526345464322643, 0.373574233602971, 0.516178154279606, 0.507524004141834 };
  vector<double> expected_resid = { 0.626153637225693, 0.219850398122767, 0.333521346184348, -0.504185638464105, 0.575670326315469, 0.157696826665612, 0.630486562300961, 0.49870923646625, 0.608591170706794, 0.329857726485879, -0.616225441901857, -0.210597307527924, -0.463290963169013, -0.697938227509708, -0.594263151190782, -0.240364041411177, 0.28249925395447, -0.385697671445548, -0.363177706176422, 0.27938680097465, 0.626382683219943, -0.0964716804290593, 0.17426391863132, -0.448680997567444, -0.521486628711504, 0.722931424597445, -0.526345464322643, -0.373574233602971, 0.483821845720394, -0.507524004141834 };
  vector<double> expected_pvals = { 0.659773503210558, 0.258673158059626, 0.0794026163233215, 0.862610307430172 };

  ASSERT_NEAR(betas[0], -0.18626, 0.00001);
  ASSERT_NEAR(betas[1], 0.55630, 0.00001);
  ASSERT_NEAR(betas[2], -0.67281, 0.00001);
  ASSERT_NEAR(betas[3], -0.07457, 0.00001);

  ASSERT_NEAR(std_err[0], 0.42311, 0.00001);
  ASSERT_NEAR(std_err[1], 0.49250, 0.00001);
  ASSERT_NEAR(std_err[2], 0.38355, 0.00001);
  ASSERT_NEAR(std_err[3], 0.43093, 0.00001);

  for (int i = 0; i < expected_p.size(); i++) {
    ASSERT_NEAR(p[i], expected_p[i], 0.00001);
  }

  for (int i = 0; i < expected_resid.size(); i++) {
    ASSERT_NEAR(resid[i], expected_resid[i], 0.00001);
  }

  for (int i = 0; i < expected_pvals.size(); i++) {
    ASSERT_NEAR(pvals[i], expected_pvals[i], 0.00001);
  }

  int foo = 1;
}

TEST_F(LDServerTest, regression_linear_small) {
  arma::mat x = {
    { 1, -0.096804537415269, 0.394521501456981, 0.136719199203261  },
    { 1, 1.87480075412442, -0.576503014867386, -0.294691393024923  },
    { 1, 0.772754744319657, -0.558625773661303, -0.976660217651494  },
    { 1, 0.717732517365932, 0.0936844055088714, 1.78662210708278  },
    { 1, 1.8018548122351, 1.85515254706929, -1.70359905461682  },
    { 1, -0.348127414820295, -3.05083075001173, -0.0369152710070338  },
    { 1, -1.01502457917055, -0.276941129747707, -0.406115851940272  },
    { 1, 1.10078194115342, 0.515343122857823, 0.995075331192894  },
    { 1, 0.34158930036517, 0.570973367742649, 0.817951168123065  },
    { 1, 0.778222729194391, -0.717535882710764, 0.276220595509114  },
    { 1, -0.445580619793771, -1.35355138636549, 0.0400984960691716  },
    { 1, 1.05968330903572, 2.48073248812764, 0.744166070387559  },
    { 1, 1.13193063729192, 0.803345354713477, 0.67083197357916  },
    { 1, 0.64438005748593, -0.820235260483927, -1.52106595803099  },
    { 1, 0.998139980346589, -0.101007365065264, 0.742100369064892  },
    { 1, 0.0043382638725237, 1.44109213606475, -0.0368974515646374  },
    { 1, 1.2213571153, -0.547334882194856, -0.947655460648014  },
    { 1, -0.545318955493577, -0.167697594444228, 1.18861767550958  },
    { 1, -0.11332905303665, 0.477938596320964, -0.124355176118656  },
    { 1, -1.06063645191269, -2.49783021515651, -0.579642166533558  },
    { 1, -1.31296784040276, -0.597392049046657, 0.0267899494002477  },
    { 1, -1.14312422332763, 2.12767104213113, -0.223249127365326  },
    { 1, 0.832580243798403, -1.78912384982458, -1.00652278049377  },
    { 1, -1.12411740497594, -0.795110424818677, -0.947332980605162  },
    { 1, 0.705237914919055, 0.222294866312317, -0.39563859497371  },
    { 1, 0.229680544111858, 1.25628142318361, 0.741735315848646  },
    { 1, -1.0101433303265, -1.23653274314961, -0.291369194678221  },
    { 1, 0.82667455981141, 1.33750612592533, -1.46654503496563  },
    { 1, 0.528798998554926, 0.190260018116963, -1.13765940889545  },
    { 1, 1.18406942883136, 0.851258176994322, -1.74868221717553  },
  };

  arma::vec y = { -0.0788360277700111, -1.07581634483181, 0.640742961590315, 0.186770707908177, -0.613528194812187, 0.960099690164341, -0.672106790288292, 0.79712388914405, 1.63303901014036, 1.11635186445792, -0.761201159517016, 1.08821794453735, -0.276746184979133, -0.391392222221473, -1.56053575864685, -1.10783780966815, -0.872155105408476, 1.47352970670923, -0.743947925727988, 1.05830336204417, 0.405522701348294, -0.141060188701722, 1.90823438037963, -0.138431608664085, 1.48514643537764, 0.10547183676043, -0.870211784368017, 1.33681171252519, -0.165743329688721, 0.232370350684276 };
  LinearRegression model;
  model.fit(y, x);

  /*
    Call:
    lm(formula = yqt ~ x)

    Residuals:
        Min      1Q  Median      3Q     Max
    -1.8923 -0.7055 -0.1251  0.6987  1.6438

    Coefficients:
                Estimate Std. Error t value Pr(>|t|)
    (Intercept)  0.17385    0.19244   0.903    0.375
    x1           0.05108    0.21299   0.240    0.812
    x2          -0.10021    0.15282  -0.656    0.518
    x3           0.13040    0.20584   0.634    0.532

    Residual standard error: 0.987 on 26 degrees of freedom
    Multiple R-squared:  0.0272,	Adjusted R-squared:  -0.08504
    F-statistic: 0.2423 on 3 and 26 DF,  p-value: 0.866
   */

  auto betas = *model.getBetas();
  auto std_err = *model.getStandardErrors();
  auto resid = *model.getResiduals();
  auto pvals = *model.getPvalues();
  double sigma2 = model.getSigmaSquared();

  vector<double> expected_betas = { 0.173849916958154, 0.0510761805307706, -0.100208889616146, 0.130403963163449 };
  vector<double> expected_se = { 0.192441744917242, 0.212990052660106, 0.152819730775099, 0.205838801429248 };
  vector<double> expected_pvals = { 0.374613288944015, 0.812361602174763, 0.517757136787951, 0.531925016699081 };
  vector<double> expected_resid = { -0.226035702524885, -1.36476572498647, 0.498804778361553, -0.247332837868399, -0.471351128276382, 0.503124327668361, -0.768905975102108, 0.488930430313695, 1.39229454959383, 0.794829568462346, -1.05315940456796, 1.0117926948211, -0.515387597513274, -0.481996446843061, -1.89226151852352, -1.13268749187078, -1.03965707189118, 1.15572735386558, -0.847899323666076, 0.763909547403557, 0.235376657362428, -0.0141995628784897, 1.64382788979687, -0.211007059794852, 1.34914442184545, -0.0509439434755223, -1.07838331365971, 1.44601170496088, -0.199181238957555, 0.311381417944568 };
  double expected_sigma2 = 0.9741761;

  for (int i = 0; i < expected_betas.size(); i++) {
    ASSERT_NEAR(betas[i], expected_betas[i], 0.00001);
  }

  for (int i = 0; i < expected_resid.size(); i++) {
    ASSERT_NEAR(resid[i], expected_resid[i], 0.00001);
  }

  for (int i = 0; i < expected_se.size(); i++) {
    ASSERT_NEAR(std_err[i], expected_se[i], 0.00001);
  }

  for (int i = 0; i < expected_pvals.size(); i++) {
    ASSERT_NEAR(pvals[i], expected_pvals[i], 0.00001);
  }

  ASSERT_NEAR(sigma2, expected_sigma2, 0.00001);

  int foo = 1;
}