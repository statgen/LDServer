#ifndef LDSERVER_SUMMARYSTATISTICSLOADER_H
#define LDSERVER_SUMMARYSTATISTICSLOADER_H

#include <string>
#include <set>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <tabixpp.hpp>
#include <boost/format.hpp>
#include "Types.h"

const uint32_t INIT_QUERY_LIMIT = 10000000;

struct ScoreColumns {
  uint16_t colChrom;
  uint16_t colPos;
  uint16_t colRef;
  uint16_t colAlt;
  uint16_t colAltFreq;
  uint16_t colU;
  uint16_t colV;
  uint16_t colEffectAllele;
  uint16_t colPvalue;
};

struct CovarianceColumns {
  uint16_t colCov;
  uint16_t colPos;
  uint16_t colChrom;
};

const ScoreColumns SCORE_COLUMNS_RVTEST {0, 1, 2, 3, 5, 12, 13, 3, 15};
const ScoreColumns SCORE_COLUMNS_RAREMETAL {0, 1, 2, 3, 5, 13, 14, 3, 16};
const CovarianceColumns COV_COLUMNS_RVTEST {5, 4, 0};
const CovarianceColumns COV_COLUMNS_RAREMETAL {3, 2, 0};

enum class ScoreCovFormat {RVTEST, RAREMETAL};

class ScoreCovParseException : public std::runtime_error { using std::runtime_error::runtime_error; };

class SummaryStatisticsLoader {
protected:
  std::string score_path;
  std::string cov_path;
  std::map<uint64_t, double> alt_freq;
  std::map<uint64_t, std::string> pos_variant;

  shared_ptr<LDQueryResult> cov_result;
  shared_ptr<ScoreStatQueryResult> score_result;

  ScoreCovFormat detected_format;
  double sigma2;
  uint64_t nsamples;

  void parseHeader(const std::string& filepath);
  static uint64_t getNumberOfVariantsFromCovFile(const std::string& filepath, const std::string& region);

  double getAltFreqForPosition(uint64_t& pos);
  std::string getVariantForPosition(uint64_t& pos);

  void load_scores(const std::string& chromosome, uint64_t start, uint64_t stop);
  void load_cov(const std::string& chromosome, uint64_t start, uint64_t stop);
public:
  SummaryStatisticsLoader(const std::string& score_path, const std::string& cov_path);
  void load_region(const std::string& chromosome, uint64_t start, uint64_t stop);
  shared_ptr<LDQueryResult> getCovResult();
  shared_ptr<ScoreStatQueryResult> getScoreResult();
  double getSigma2() { return sigma2; }
  uint64_t getNumSamples() { return nsamples; }
};

#endif //LDSERVER_SUMMARYSTATISTICSLOADER_H
