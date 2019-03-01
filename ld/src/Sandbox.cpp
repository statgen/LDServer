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
#include <armadillo>
using namespace std;

void test1() {
  LDServer ld_server(100);
  LDQueryResult ld_result(1000);

  ScoreServer score_server(100);
  ScoreStatQueryResult score_results(1000);
  auto segments = make_shared_segment_vector();

  string genotype_file = "../../../data/chr22.test.vcf.gz";
  ld_server.set_file(genotype_file);

  score_server.set_genotypes_file(genotype_file, 1);

  ColumnTypeMap ctmap;
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  string phenotype_file = "../../../data/chr22.test.tab";
  score_server.load_phenotypes_file(phenotype_file, ctmap, 2504, "\t", "iid", 1);
  score_server.set_phenotype("rand_qt");

  string chrom = "22";
  auto start = 50276998ul;
  auto stop = 50357719ul;

  // try out mask
  Mask mask("../../../data/mask.epacts.chr22.gencode-exons-AF01.tab.gz", "AF < 0.01", VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, chrom, start, stop);

  // try out runner
  vector<Mask> masks;
  masks.emplace_back(mask);

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

  config->pprint();

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

int main() {
  test2();
  return 0;
}