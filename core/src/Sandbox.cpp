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

int main() {
  test1();
  //test3();
  return 0;
}