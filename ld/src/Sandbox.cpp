// This program is entirely for quick testing/debugging and not meant to be used in production.

#include <string>
#include <vector>
#include "Types.h"
#include "LDServer.h"
#include "ScoreServer.h"
#include "Phenotypes.h"
#include "Mask.h"
#include "ScoreCovarianceRunner.h"
#include <armadillo>
using namespace std;

int main() {
  LDServer ld_server(100);
  LDQueryResult ld_result(1000);

  ScoreServer score_server(100);
  ScoreStatQueryResult score_results(1000);
  auto segments = make_shared_segment_vector();

  string genotype_file = "../../../data/chr22.test.vcf.gz";
  ld_server.set_file(genotype_file);

  score_server.set_genotypes_file(genotype_file, 1);

  ColumnTypeMap ctmap;
  ctmap.add("fid", ColumnType::TEXT);
  ctmap.add("iid", ColumnType::TEXT);
  ctmap.add("patid", ColumnType::TEXT);
  ctmap.add("matid", ColumnType::TEXT);
  ctmap.add("sex", ColumnType::CATEGORICAL);
  ctmap.add("rand_binary", ColumnType::CATEGORICAL);
  ctmap.add("rand_qt", ColumnType::FLOAT);

  string phenotype_file = "../../../data/chr22.test.tab";
  score_server.load_phenotypes_file(phenotype_file, ctmap, 2504, 1);
  score_server.set_phenotype("rand_qt");

  // try out mask
  Mask mask("../../../data/mask.epacts.chr22.gencode-exons-AF01.tab.gz", "AF < 0.01", VariantGroupType::GENE, GroupIdentifierType::ENSEMBL, "22", 50276998ul, 50357719ul);

  // try out runner
  vector<Mask> masks;
  masks.emplace_back(mask);

  auto config = make_score_covariance_config();

  config->segment_size = 1000;
  config->chrom = "22";
  config->start = 50276998ul;
  config->stop = 50357719ul;
  config->masks = masks;
  config->sample_subset = "ALL";
  config->genotype_files = {genotype_file};
  config->genotype_dataset_id = 1;
  config->phenotype_file = phenotype_file;
  config->column_types = ctmap;
  config->phenotype_dataset_id = 1;
  config->phenotype = "rand_qt";
  config->nrows = 2504;

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