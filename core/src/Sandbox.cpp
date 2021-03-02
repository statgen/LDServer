// This program is entirely for quick testing/debugging and not meant to be used in production.

#include <string>
#include <vector>
#include "Types.h"
#include "LDServer.h"
#include "ScoreServer.h"
#include "Phenotypes.h"
#include "Mask.h"
#include "ScoreCovarianceRunner.h"
#include "Raw.h"
#include <chrono>
#include <armadillo>
using namespace std;

void test1() {
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
  config->genotype_files = {"../../../data/chr22.monomorphic_test.vcf.gz"};
  config->genotype_dataset_id = 1;
  config->phenotype_file = "../../../data/chr22.test.missing_values.tab";
  config->phenotype_column_types = ctmap;
  config->phenotype_dataset_id = 1;
  config->phenotype = "rand_qt";
  config->phenotype_nrows = 2504;
  config->phenotype_sample_column = "iid";
  config->phenotype_delim = "\t";
  config->pprint();

  // Load mask
  Mask mask("../../../data/mask.epacts.chr22.gencode-exons-AF01.tab.gz", 1, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, chrom, start, stop);
  vector<Mask> masks;
  masks.emplace_back(mask);
  config->masks = masks;

  // Execute runner
  ScoreCovarianceRunner runner(config);
  runner.run();

  string json = runner.getPrettyJSON();
  //cout << json << endl;

  ofstream out("json.txt");
  out << json;
  out.close();

  int x = 0;
}

void test2() {
  auto samples1 = extract_samples("../../../data/chr21.test.vcf.gz");
  auto samples2 = extract_samples("../../../data/chr21.test.sav");
}

void test3() {
  auto time_start = std::chrono::system_clock::now();

  // Genotype spec
  uint32_t genotype_dataset_id = 2;
  string genotype_file = "../../../data/METSIM/data/genotypes_chr22.vcf.gz";

  // Phenotype spec
  ColumnTypeMap ctmap;
  ctmap.add("studyid", ColumnType::TEXT);
  ctmap.add("hba1c", ColumnType::FLOAT);
  ctmap.add("ldl", ColumnType::FLOAT);
  ctmap.add("HOMA", ColumnType::FLOAT);
  ctmap.add("t2d", ColumnType::CATEGORICAL);
  ctmap.add("hdl", ColumnType::FLOAT);
  ctmap.add("sbp", ColumnType::FLOAT);
  ctmap.add("ins_30m", ColumnType::FLOAT);
  ctmap.add("ins_fast", ColumnType::FLOAT);
  ctmap.add("HOMA_B", ColumnType::FLOAT);
  ctmap.add("coffee", ColumnType::FLOAT);
  ctmap.add("glu_2h", ColumnType::FLOAT);
  ctmap.add("chol", ColumnType::FLOAT);
  ctmap.add("matsuda", ColumnType::FLOAT);
  ctmap.add("glu_fast", ColumnType::FLOAT);
  ctmap.add("ins_2h", ColumnType::FLOAT);
  ctmap.add("tg", ColumnType::FLOAT);
  ctmap.add("glu_30m", ColumnType::FLOAT);
  ctmap.add("dbp", ColumnType::FLOAT);

  string phenotype = "ldl";
  uint64_t phenotype_nrows = 9376;
  string phenotype_sample_col = "studyid";
  uint32_t phenotype_dataset_id = 3;
  string phenotype_file = "../../../data/METSIM/data/metsim.for_raremetal.pheno.tab";

  // Region
  string chrom = "22";
  auto start = 50276998ul;
  auto stop = 50357719ul;

  // Masks for runner
  uint64_t mask_id = 3;
  Mask mask("../../../data/METSIM/data/mask.epacts.allchroms.gencode-exons-AF01.tab.gz", mask_id, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, chrom, start, stop);
  vector<Mask> masks;
  masks.emplace_back(mask);

  // Config for runner
  auto config = make_score_covariance_config();
  config->chrom = chrom;
  config->start = start;
  config->stop = stop;
  config->segment_size = 1000;
  config->masks = masks;
  config->sample_subset = "ALL";
  config->genotype_files = {genotype_file};
  config->genotype_dataset_id = genotype_dataset_id;
  config->phenotype_file = phenotype_file;
  config->phenotype_column_types = ctmap;
  config->phenotype_dataset_id = phenotype_dataset_id;
  config->phenotype = phenotype;
  config->phenotype_nrows = phenotype_nrows;
  config->phenotype_sample_column = phenotype_sample_col;
  config->phenotype_delim = "\t";

  config->pprint();

  ScoreCovarianceRunner runner(config);
  runner.run();

  string json = runner.getPrettyJSON();
  //cout << json << endl;

  ofstream out("json.txt");
  out << json;
  out.close();

  std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - time_start;
  cout << "Time required: " << elapsed.count() << endl;
}

void test4() {
  auto time_start = std::chrono::system_clock::now();

  // Genotype spec
  uint32_t genotype_dataset_id = 2;
  string genotype_file = "../../../data/METSIM/data/genotypes_chr22.vcf.gz";

  // Phenotype spec
  ColumnTypeMap ctmap;
  ctmap.add("studyid", ColumnType::TEXT);
  ctmap.add("hba1c", ColumnType::FLOAT);
  ctmap.add("ldl", ColumnType::FLOAT);
  ctmap.add("HOMA", ColumnType::FLOAT);
  ctmap.add("t2d", ColumnType::CATEGORICAL);
  ctmap.add("hdl", ColumnType::FLOAT);
  ctmap.add("sbp", ColumnType::FLOAT);
  ctmap.add("ins_30m", ColumnType::FLOAT);
  ctmap.add("ins_fast", ColumnType::FLOAT);
  ctmap.add("HOMA_B", ColumnType::FLOAT);
  ctmap.add("coffee", ColumnType::FLOAT);
  ctmap.add("glu_2h", ColumnType::FLOAT);
  ctmap.add("chol", ColumnType::FLOAT);
  ctmap.add("matsuda", ColumnType::FLOAT);
  ctmap.add("glu_fast", ColumnType::FLOAT);
  ctmap.add("ins_2h", ColumnType::FLOAT);
  ctmap.add("tg", ColumnType::FLOAT);
  ctmap.add("glu_30m", ColumnType::FLOAT);
  ctmap.add("dbp", ColumnType::FLOAT);

  string phenotype = "ldl";
  uint64_t phenotype_nrows = 9376;
  string phenotype_sample_col = "studyid";
  uint32_t phenotype_dataset_id = 3;
  string phenotype_file = "../../../data/METSIM/data/metsim.for_raremetal.pheno.tab";

  // Region
  string chrom = "22";
  auto start = 50276998ul;
  auto stop = 50357719ul;

  // Masks for runner
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

  // Config for runner
  auto config = make_score_covariance_config();
  config->chrom = chrom;
  config->start = start;
  config->stop = stop;
  config->segment_size = 1000;
  config->masks = masks;
  config->sample_subset = "ALL";
  config->genotype_files = {genotype_file};
  config->genotype_dataset_id = genotype_dataset_id;
  config->phenotype_file = phenotype_file;
  config->phenotype_column_types = ctmap;
  config->phenotype_dataset_id = phenotype_dataset_id;
  config->phenotype = phenotype;
  config->phenotype_nrows = phenotype_nrows;
  config->phenotype_sample_column = phenotype_sample_col;
  config->phenotype_delim = "\t";

  config->pprint();

  ScoreCovarianceRunner runner(config);
  runner.run();

  string json = runner.getPrettyJSON();
  //cout << json << endl;

  ofstream out("json.txt");
  out << json;
  out.close();

  std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - time_start;
  cout << "Time required: " << elapsed.count() << endl;
}

void perf_sav_55k() {
  auto time_start = std::chrono::system_clock::now();

  string genotype_file = "../../../private/55k-exomes/55k.clean.chrX.sav";
  string phenotype_file = "../../../private/55k-exomes/55kQTsRemoved.ped";

  ColumnTypeMap ctmap;
  ctmap.add("fid", ColumnType::TEXT);
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("patid", ColumnType::TEXT);
  ctmap.add("matid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("t2d", ColumnType::CATEGORICAL);
  ctmap.add("Age", ColumnType::FLOAT);
  ctmap.add("DBP_ADJ", ColumnType::FLOAT);
  ctmap.add("DBP", ColumnType::FLOAT);
  ctmap.add("SBP_ADJ", ColumnType::FLOAT);
  ctmap.add("SBP", ColumnType::FLOAT);
  ctmap.add("FAST_GLU_ADJ", ColumnType::FLOAT);
  ctmap.add("FAST_GLU", ColumnType::FLOAT);
  ctmap.add("FAST_INS_ADJ", ColumnType::FLOAT);
  ctmap.add("FAST_INS", ColumnType::FLOAT);
  ctmap.add("HDL_ADJ", ColumnType::FLOAT);
  ctmap.add("HDL", ColumnType::FLOAT);
  ctmap.add("LDL_ADJ", ColumnType::FLOAT);
  ctmap.add("LDL", ColumnType::FLOAT);
  ctmap.add("BMI_ADJ", ColumnType::FLOAT);
  ctmap.add("BMI", ColumnType::FLOAT);
  ctmap.add("SERUM_CREATININE_ADJ", ColumnType::FLOAT);
  ctmap.add("SERUM_CREATININE", ColumnType::FLOAT);
  ctmap.add("HR2_GLU_ADJ", ColumnType::FLOAT);
  ctmap.add("HR2_GLU", ColumnType::FLOAT);
  ctmap.add("HR2_INS_ADJ", ColumnType::FLOAT);
  ctmap.add("HR2_INS", ColumnType::FLOAT);
  ctmap.add("CHOL_ADJ", ColumnType::FLOAT);
  ctmap.add("CHOL", ColumnType::FLOAT);
  ctmap.add("TG_ADJ", ColumnType::FLOAT);
  ctmap.add("TG", ColumnType::FLOAT);
  ctmap.add("ANCESTRY", ColumnType::CATEGORICAL);
  ctmap.add("origin", ColumnType::CATEGORICAL);
  ctmap.add("C1", ColumnType::FLOAT);
  ctmap.add("C2", ColumnType::FLOAT);
  ctmap.add("C3", ColumnType::FLOAT);
  ctmap.add("C4", ColumnType::FLOAT);
  ctmap.add("C5", ColumnType::FLOAT);
  ctmap.add("C6", ColumnType::FLOAT);
  ctmap.add("C7", ColumnType::FLOAT);
  ctmap.add("C8", ColumnType::FLOAT);
  ctmap.add("C9", ColumnType::FLOAT);
  ctmap.add("C10", ColumnType::FLOAT);

  // Region to analyze
  string chrom = "8";
  auto start = 31144764ul;
  auto stop = 33146263ul;

  // Setup mask (this would be user defined client-side via API)
  uint64_t mask_id = 0;
  vector<VariantGroup> vg_vec;
  VariantGroup vg;
  vg.chrom = "X";
  vg.name = "TEST";
  vg.start = 31144764;
  vg.stop = 33146263;
  vg.variants = {
    VariantMeta("X:33146263_C/A"),
    VariantMeta("X:31196048_C/T"),
    VariantMeta("X:32662370_T/TTC"),
    VariantMeta("X:32459433_T/A"),
    VariantMeta("X:31144764_T/A"),
    VariantMeta("X:32662247_A/C"),
    VariantMeta("X:32429867_G/A"),
    VariantMeta("X:32591646_C/T")
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
  config->phenotype = "FAST_GLU";
  config->phenotype_nrows = 43125;
  config->phenotype_sample_column = "iid";
  config->phenotype_delim = "\t";

  // Run score/covariance calculations
  string json;
  for (int i = 0; i < 10; i++) {
    ScoreCovarianceRunner runner(config);
    runner.run();
    json = runner.getJSON();
  }

  std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - time_start;
  cout << "Time required: " << elapsed.count() << endl;

  // Parse back out JSON
  rapidjson::Document doc;
  doc.Parse(json.c_str());

  int x = 0;
  // Tests
//  ASSERT_EQ(doc["data"]["groups"][0]["variants"].Size(), 18);
//  ASSERT_EQ(doc["data"]["groups"][0]["covariance"].Size(), 171);
//  ASSERT_EQ(doc["data"]["variants"][0]["variant"], "22:50354416_G/C");
//  ASSERT_NEAR(doc["data"]["variants"][0]["score"].GetDouble(), -45.31790009565313, 0.0001);
//  ASSERT_NEAR(doc["data"]["groups"][0]["covariance"][0].GetDouble(), 0.39843530436, 0.0001);
}

void mask_segfault() {
  string mask_path = "/Users/welchr/projects/LDServer/rest/raremetal/../../data/test.smallchunk.mask.epacts.tab";
  auto itype = GroupIdentifierType::ENSEMBL;
  auto gtype = VariantGroupType::GENE;
  Mask mask(
    mask_path,
    3,
    gtype,
    itype,
    "1",
    2,
    307
  );
}

void perf_ld_server() {
  auto time_start = std::chrono::system_clock::now();

  LDServer server(100);
  LDQueryResult result(100000);
  server.set_file("../../../data/chr22.test.sav");

  string json;
  for (int i = 0; i < 30; i++) {
    result.erase();
    server.compute_region_ld("22", 50244251, 51244237, correlation::LD_RSQUARE, result);

    // Get JSON
    json = result.get_json("blah");
  }

  // Parse back out JSON
  rapidjson::Document doc;
  doc.Parse(json.c_str());

  std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - time_start;
  cout << "Time required: " << elapsed.count() << endl;
}

void sumstats() {
//  SummaryStatisticsLoader loader(
//    "../../../data/test.smallchunk.MetaScore.assoc.gz",
//    "../../../data/test.smallchunk.MetaCov.assoc.gz"
//  );
//  loader.load_region("1", 2, 9);

  string chrom = "1";
  auto start = 2ul;
  auto stop = 307ul;

  uint64_t mask_id = 3;
  Mask mask("../../../data/test.smallchunk.mask.epacts.tab.gz", mask_id, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, chrom, start, stop);
  vector<Mask> masks;
  masks.emplace_back(mask);

  auto config = make_score_covariance_config();
  config->chrom = chrom;
  config->start = start;
  config->stop = stop;
  config->segment_size = 10;
  config->masks = masks;
  config->summary_stat_dataset_id = 1;
  config->summary_stat_score_files = {"../../../data/test.smallchunk.MetaScore.assoc.gz"};
  config->summary_stat_cov_files = {"../../../data/test.smallchunk.MetaCov.assoc.gz"};

  // Run score/covariance calculations
  ScoreCovarianceRunner runner(config);
  runner.run();
  auto json = runner.getJSON();
  rapidjson::Document doc;
  doc.Parse(json.c_str());
  auto& data = doc["data"];
  auto& groups = doc["data"]["groups"];
  auto& variants = doc["data"]["variants"];
  auto& group1 = doc["data"]["groups"][0];
  auto& group2 = doc["data"]["groups"][1];
  auto& group1_variant1 = doc["data"]["groups"][0]["variants"][0];
  auto& group2_variant1 = doc["data"]["groups"][1]["variants"][0];
  cout << "DONE!" << endl;
}

void tabixpp_error() {
  string chrom = "22";
  auto start = 1ul;
  auto stop = 2ul;

  uint64_t mask_id = 0;
  vector<VariantGroup> vg_vec;
  VariantGroup vg;
  vg.chrom = "22";
  vg.name = "PIM3";
  vg.start = 1;
  vg.stop = 1;
  vg.variants = {VariantMeta("1:1_A/T")};
  vg_vec.push_back(vg);
  Mask mask(mask_id, VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, vg_vec);
  vector<Mask> masks;
  masks.emplace_back(mask);

  auto config = make_score_covariance_config();
  config->chrom = chrom;
  config->start = start;
  config->stop = stop;
  config->segment_size = 1000;
  config->masks = masks;
  config->summary_stat_dataset_id = 1;
  config->summary_stat_score_files = {"/Users/welchr/projects/LDServer/rest/raremetal/../../data/test_sumstat_loader_rm.scores.assoc.gz"};
  config->summary_stat_cov_files = {"/Users/welchr/projects/LDServer/rest/raremetal/../../data/test_sumstat_loader_rm.cov.assoc.gz"};

  // Run score/covariance calculations
  ScoreCovarianceRunner runner(config);
  runner.run();
//  auto json = runner.getJSON();
//  rapidjson::Document doc;
//  doc.Parse(json.c_str());
//  auto& data = doc["data"];
//  auto& groups = doc["data"]["groups"];
//  auto& variants = doc["data"]["variants"];
//  auto& group1 = doc["data"]["groups"][0];
//  auto& group2 = doc["data"]["groups"][1];
//  auto& group1_variant1 = doc["data"]["groups"][0]["variants"][0];
//  auto& group2_variant1 = doc["data"]["groups"][1]["variants"][0];
  cout << "DONE!" << endl;
}

int main() {
  //perf_sav_55k();
  //sumstats();
  //tabixpp_error();
  //test3();
  perf_sav_55k();
  return 0;
}