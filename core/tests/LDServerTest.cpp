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
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "../src/LDServer.h"
#include "../src/ScoreServer.h"
#include "RareMetal.h"
#include "RvTest.h"
#include "../src/Phenotypes.h"
#include "../src/Mask.h"
#include "../src/ScoreCovarianceRunner.h"
#include "../src/SummaryStatisticsLoader.h"
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

    auto load_rvtest_scores(const string &file) {
      return make_shared<RvTestScores>(file);
    }

    void load_raremetal_covariance(const string &path, map<string, double> &values) {
        unique_ptr<istream> file;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        ifstream fs(path, ios_base::in | ios_base::binary);

        bool is_gz = path.find(".gz") != string::npos;
        string line;
        if (is_gz) {
          inbuf.push(boost::iostreams::gzip_decompressor());
          inbuf.push(fs);
          file = make_unique<istream>(&inbuf);
        }
        else {
          file = make_unique<ifstream>(path);
        }

        vector<string> tokens;
        vector<uint64_t> positions;
        vector<double> cov;
        auto line_separator = regex("[ \t]");
        auto field_separator = regex(",");
        string comment = "#";

        getline(*file, line); // skip header;
        while (getline(*file, line)) {
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

  void load_rvtest_covariance(const string &path, map<string, double> &values) {
    unique_ptr<istream> file;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    ifstream fs(path, ios_base::in | ios_base::binary);

    bool is_gz = path.find(".gz") != string::npos;
    string line;
    if (is_gz) {
      inbuf.push(boost::iostreams::gzip_decompressor());
      inbuf.push(fs);
      file = make_unique<istream>(&inbuf);
    }
    else {
      file = make_unique<ifstream>(path);
    }

    vector<string> tokens;
    vector<uint64_t> positions;
    vector<double> cov;
    auto line_separator = regex("[ \t]");
    auto field_separator = regex(",");
    string comment = "#";
    auto regex_header = regex("CHROM\tSTART_POS.*");

    getline(*file, line); // skip header;
    while (getline(*file, line)) {
      smatch match;
      if (std::equal(comment.begin(), comment.end(), line.begin())) {
        continue;
      }
      else if (regex_search(line, match, regex_header)) {
        continue;
      }

      // Split line and insert into tokens
      copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

      string chrom(tokens.at(0));
      string ref_pos(tokens.at(1));

      // Load the positions on this row
      transform(
        sregex_token_iterator(tokens.at(4).begin(), tokens.at(4).end(), field_separator, -1),
        sregex_token_iterator(),
        back_inserter(positions),
        [](const string &str) { return stoi(str); }
      );

      // Load the covariance values on this row
      transform(
        sregex_token_iterator(tokens.at(5).begin(), tokens.at(5).end(), field_separator, -1),
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

TEST_F(LDServerTest, summary_stat_load_raremetal_test) {
  // Load from disk using our new summary stat loader
  SummaryStatisticsLoader loader({"test_sumstat_loader_rm.scores.assoc.gz"}, {"test_sumstat_loader_rm.cov.assoc.gz"});
  loader.load_region("22", 27021502, 27026606);

  // Use our testing methods to load the same data for later comparison
  map<string, double> gold_cov;
  this->load_raremetal_covariance("test_sumstat_loader_rm.cov.assoc.gz", gold_cov);
  auto gold_scores = load_raremetal_scores("test_sumstat_loader_rm.scores.assoc.gz");

  // Check covariance
  auto cov_result = loader.getCovResult();
  ASSERT_FALSE(cov_result->data.empty());
  for (auto&& entry : cov_result->data) {
    string key(to_string(entry.position1) + "_" + to_string(entry.position2));
    double value_gold = abs(gold_cov.find(key)->second); // needs abs because the loader flips cov according to alt allele freq
    double value_test = abs(entry.value);
    ASSERT_NE(entry.variant1, "");
    ASSERT_NE(entry.variant2, "");
    ASSERT_NEAR(value_gold, value_test, 0.0001);
  }

  // Check scores
  auto score_result = loader.getScoreResult();
  ASSERT_FALSE(score_result->data.empty());
  for (auto&& score_res : score_result->data) {
    auto gold = gold_scores->get_record(score_res.variant);
    ASSERT_NEAR(gold->u_stat, score_res.score_stat, 0.001);
    ASSERT_NEAR(gold->pvalue, score_res.pvalue, 0.001);
  }
  ASSERT_TRUE(loader.getNumSamples() > 0);
}

TEST_F(LDServerTest, summary_stat_load_rvtest_test) {
  // Load from disk using our new summary stat loader
  SummaryStatisticsLoader loader({"test.smallchunk.MetaScore.assoc.gz"}, {"test.smallchunk.MetaCov.assoc.gz"});
  loader.load_region("1", 2, 307);

  // Use our testing methods to load the same data for later comparison
  map<string, double> gold_cov;
  this->load_rvtest_covariance("test.smallchunk.MetaCov.assoc.gz", gold_cov);
  auto gold_scores = load_rvtest_scores("test.smallchunk.MetaScore.assoc.gz");

  // Check covariance
  auto cov_result = loader.getCovResult();
  ASSERT_FALSE(cov_result->data.empty());
  for (auto&& entry : cov_result->data) {
    string key(to_string(entry.position1) + "_" + to_string(entry.position2));
    double value_gold = gold_cov.find(key)->second;
    double value_test = entry.value;
    ASSERT_NE(entry.variant1, "");
    ASSERT_NE(entry.variant2, "");
    ASSERT_NEAR(value_gold, value_test, 0.0001);
  }

  // Check scores
  auto score_result = loader.getScoreResult();
  ASSERT_FALSE(score_result->data.empty());
  for (auto&& score_res : score_result->data) {
    auto gold = gold_scores->get_record(score_res.variant);
    ASSERT_NEAR(gold->u_stat, score_res.score_stat, 0.001);
    ASSERT_NEAR(gold->pvalue, score_res.pvalue, 0.001);
  }
  ASSERT_TRUE(loader.getNumSamples() > 0);
}

TEST_F(LDServerTest, summary_stat_load_rvtest_bad_ustat) {
  using ::testing::HasSubstr;

  // Load from disk using our new summary stat loader
  SummaryStatisticsLoader loader({"rvtest_score_fail_ustat.gz"}, {"rvtest_cov_fail_base.gz"});

  try {
    loader.load_region("9", 22133, 22133);
    FAIL() << "Expected LDServerGenericException";
  }
  catch(LDServerGenericException& e) {
    EXPECT_THAT(e.get_secret(), HasSubstr("offending value was 'NA' in column 'U_STAT' for variant '9:22133_A/T'"));
  }
  catch(...) {
    FAIL() << "Expected LDServerGenericException";
  }
}

TEST_F(LDServerTest, summary_stat_load_rvtest_missing_af) {
  // Load from disk using our new summary stat loader
  SummaryStatisticsLoader loader({"test.afmissing.MetaScore.assoc.gz"}, {"test.smallchunk.MetaCov.assoc.gz"});
  loader.load_region("1", 2, 2);

  auto score_res = loader.getScoreResult();
  ASSERT_NEAR(score_res->data[0].alt_freq, 0.004464285, 0.00001);
}

TEST_F(LDServerTest, summary_stat_load_rvtest_noheader_test) {
  // Load from disk using our new summary stat loader
  SummaryStatisticsLoader loader({"test.smallchunk.noheader.MetaScore.assoc.gz"}, {"test.smallchunk.noheader.MetaCov.assoc.gz"});
  loader.load_region("1", 2, 307);

  // Use our testing methods to load the same data for later comparison
  map<string, double> gold_cov;
  this->load_rvtest_covariance("test.smallchunk.MetaCov.assoc.gz", gold_cov);
  auto gold_scores = load_rvtest_scores("test.smallchunk.MetaScore.assoc.gz");

  // Check covariance
  auto cov_result = loader.getCovResult();
  ASSERT_FALSE(cov_result->data.empty());
  for (auto&& entry : cov_result->data) {
    string key(to_string(entry.position1) + "_" + to_string(entry.position2));
    double value_gold = gold_cov.find(key)->second;
    double value_test = entry.value;
    ASSERT_NE(entry.variant1, "");
    ASSERT_NE(entry.variant2, "");
    ASSERT_NEAR(value_gold, value_test, 0.0001);
  }

  // Check scores
  auto score_result = loader.getScoreResult();
  ASSERT_FALSE(score_result->data.empty());
  for (auto&& score_res : score_result->data) {
    auto gold = gold_scores->get_record(score_res.variant);
    ASSERT_NEAR(gold->u_stat, score_res.score_stat, 0.001);
    ASSERT_NEAR(gold->pvalue, score_res.pvalue, 0.001);
  }
  ASSERT_TRUE(loader.getNumSamples() > 0);
}

TEST_F(LDServerTest, summary_stat_ldserver_compare) {
  // Summary stat loader
  SummaryStatisticsLoader loader({"chr21.test.RAND_QT.singlevar.score.txt.gz"}, {"chr21.test.RAND_QT.singlevar.cov.txt.gz"});
  loader.load_region("21", 9411239, 9411793);

  auto loader_scores = loader.getScoreResult();
  auto loader_cov = loader.getCovResult();
  auto loader_sigma2 = loader.getSigma2();
  auto loader_samples = loader.getNumSamples();

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
  ASSERT_NEAR(loader_sigma2, score_result.sigma2, 0.001);
  ASSERT_TRUE(loader_samples == score_result.nsamples);

  ASSERT_TRUE(loader_scores->data.size() == score_result.data.size());
  ASSERT_TRUE(loader_cov->data.size() == ld_result.data.size());

  for (int i = 0; i < loader_scores->data.size(); i++) {
    double u_loader = loader_scores->data[i].score_stat;
    double u_server = score_result.data[i].score_stat;
    ASSERT_NEAR(u_loader, u_server, 0.001);
  }

  loader_cov->sort_by_variant();
  ld_result.sort_by_variant();
  for (int i = 0; i < loader_cov->data.size(); i++) {
    auto& result_loader = loader_cov->data[i];
    auto& result_server = ld_result.data[i];

    ASSERT_TRUE(result_loader.variant1 == result_server.variant1);
    ASSERT_TRUE(result_loader.variant2 == result_server.variant2);

    double v_loader = result_loader.value * loader_sigma2;
    double v_server = result_server.value;
    ASSERT_NEAR(v_loader, v_server, 0.001);
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
    FAIL() << "Expected LDServerGenericException";
  }
  catch (LDServerGenericException& e) {
    EXPECT_THAT(e.what(), HasSubstr("Error reading phenotype file on line 3, column 3 (rand_qt)"));
    EXPECT_THAT(e.what(), HasSubstr("invalid value: something"));
  }
  catch (...) {
    FAIL() << "Expected LDServerGenericException";
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

TEST_F(LDServerTest, pheno_no_exist) {
  Phenotypes pheno;
  ColumnTypeMap ctmap;
  ASSERT_THROW(pheno.load_file("does-not-exist.ped", ctmap, 2504, "\t", "iid"), LDServerGenericException);
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
