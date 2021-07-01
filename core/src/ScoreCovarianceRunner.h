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
#include "SummaryStatisticsLoader.h"
#include "RaremetalSummaryStatisticsLoader.h"
#include "MetastaarSummaryStatisticsLoader.h"
#include "Phenotypes.h"
#include "Segment.h"
#include "Types.h"

enum class VariantFormat {EPACTS, COLONS};

class ScoreCovarianceConfig {
public:
  /**
   * Region specification
   */
  std::string chrom;
  uint64_t start;
  uint64_t stop;

  /**
   * Relevant settings when genotype and phenotype files are specified.
   */
  std::vector<std::string> genotype_files;
  uint32_t genotype_dataset_id;
  std::string phenotype_file;
  uint32_t phenotype_dataset_id;
  std::string phenotype;
  ColumnTypeMap phenotype_column_types;
  std::vector<std::string> phenotype_analysis_columns;
  uint64_t phenotype_nrows;
  std::string phenotype_delim;
  std::string phenotype_sample_column;
  std::string sample_subset;
  std::vector<std::string> samples;

  /**
   * Settings for when serving scores/covariance from rvtest or raremetalworker generated files
   */
  uint32_t summary_stat_dataset_id;
  std::vector<std::string> summary_stat_score_files;
  std::vector<std::string> summary_stat_cov_files;
  std::string summary_stat_format;

  /**
   * Mask related settings
   */
  std::vector<Mask> masks;

  /**
   * Cache related settings
   */
  uint32_t segment_size;
  std::string redis_hostname;
  uint16_t redis_port;

  /**
   * Output related settings
   */
  VariantFormat variant_format = VariantFormat::EPACTS;

  void pprint() const;
};

shared_ptr<ScoreCovarianceConfig> make_score_covariance_config();

enum class ScoreCovRunMode {COMPUTE, PRECOMPUTE};

class ScoreCovarianceRunner {
protected:
  std::shared_ptr<rapidjson::Document> document;
  std::shared_ptr<ScoreCovarianceConfig> config;
  std::shared_ptr<LDServer> ld_server;
  std::shared_ptr<ScoreServer> score_server;
  std::shared_ptr<SummaryStatisticsLoader> summary_stat_loader;
  ScoreCovRunMode run_mode;
public:
  ScoreCovarianceRunner(std::shared_ptr<ScoreCovarianceConfig> config);
  void run();
  std::string getJSON() const;
  std::string getPrettyJSON() const;
};

#endif //LDSERVER_RAREMETALRUNNER_H
