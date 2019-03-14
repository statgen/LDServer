#ifndef LDSERVER_RAREMETALRUNNER_H
#define LDSERVER_RAREMETALRUNNER_H

#include <iostream>
#include <boost/format.hpp>
#include <memory>
#include <string>
#include <vector>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/prettywriter.h>
#include <cereal/external/rapidjson/writer.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include "Mask.h"
#include "LDServer.h"
#include "ScoreServer.h"
#include "Phenotypes.h"
#include "Segment.h"
#include "Types.h"

class ScoreCovarianceConfig {
public:
  std::string chrom;
  uint64_t start;
  uint64_t stop;
  std::vector<std::string> genotype_files;
  uint32_t genotype_dataset_id;
  std::string phenotype_file;
  uint32_t phenotype_dataset_id;
  std::string phenotype;
  ColumnTypeMap phenotype_column_types;
  uint64_t phenotype_nrows;
  std::string phenotype_delim;
  std::string phenotype_sample_column;
  std::vector<Mask> masks;
  std::string sample_subset;
  std::vector<std::string> samples;
  uint32_t segment_size;
  std::string redis_hostname;
  uint16_t redis_port;

  void pprint() const;
};

shared_ptr<ScoreCovarianceConfig> make_score_covariance_config();

class ScoreCovarianceRunner {
protected:
  std::shared_ptr<rapidjson::Document> document;
  std::shared_ptr<ScoreCovarianceConfig> config;
public:
  ScoreCovarianceRunner(std::shared_ptr<ScoreCovarianceConfig> config);
  void run();
  std::string getJSON() const;
  std::string getPrettyJSON() const;
};

#endif //LDSERVER_RAREMETALRUNNER_H
