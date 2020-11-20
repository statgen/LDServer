#ifndef LDSERVER_SUMMARYSTATISTICSLOADER_H
#define LDSERVER_SUMMARYSTATISTICSLOADER_H

#include <string>
#include <set>
#include <vector>
#include <map>
#include <memory>
#include <limits>
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

class NoVariantsInRange : public std::runtime_error { using std::runtime_error::runtime_error; };

/**
 * Loader for "summary statistic" datasets. These are comprised of:
 *
 *   1. A file of score statistics per single variant (along with other information such as p-value, effect size, etc.)
 *   2. A file of the covariance of the score statistics for each pair of variants (usually in a sliding window across the genome)
 *
 * These files are usually generated by either rvtests, or RAREMETALWORKER. They are normally used for conducting either
 * single variant meta-analysis, or meta-analysis of gene-based aggregation tests such as the burden or SKAT test.
 *
 * https://github.com/zhanxw/rvtests#meta-analysis-models
 * https://genome.sph.umich.edu/wiki/RAREMETALWORKER
 */
class SummaryStatisticsLoader {
protected:
  std::map<std::string, std::string> score_map;
  std::map<std::string, std::string> cov_map;
  std::map<uint64_t, double> alt_freq;
  std::map<uint64_t, std::string> pos_variant;

  shared_ptr<LDQueryResult> cov_result;
  shared_ptr<ScoreStatQueryResult> score_result;

  ScoreCovFormat detected_format;
  double sigma2 = numeric_limits<double>::quiet_NaN();
  uint64_t nsamples = numeric_limits<double>::quiet_NaN();

  /**
   * Function that parses the score statistic file's header to extract:
   *    1. program name - the program (rvtest, raremetalworker) was used to create these statistics, which ends up
   *       stored in a member variable `detected_format`
   *    2. sigma2 - residual variance under the null model
   *    3. nsamples - number of samples analyzed in the model
   * @param filepath
   */
  void parseHeader(const std::string& filepath);

  /**
   * Currently unused function for getting number of variants in a covariance file, for the purposes of allocating
   * a matrix ahead of time. We currently do not need to form the matrix unless future features are required.
   * @param filepath Path to covariance matrix file.
   * @param region Region of the file to extract, given as chrom:start-end.
   * @return Number of variants within the given region in the file.
   */
  static uint64_t getNumberOfVariantsFromCovFile(const std::string& filepath, const std::string& region);

  double getAltFreqForPosition(uint64_t& pos);
  std::string getVariantForPosition(uint64_t& pos);

  /**
   * Load score statistics from a file in a given region. Same idea with load_cov().
   * These functions are both called automatically by load_region().
   * @param chromosome
   * @param start
   * @param stop
   */
  void load_scores(const std::string& chromosome, uint64_t start, uint64_t stop);
  void load_cov(const std::string& chromosome, uint64_t start, uint64_t stop);
public:
  /**
   * Create a loader object.
   * @param score_path Path on disk to the score statistic file. Should be bgzipped and tabix-indexed.
   * @param cov_path Path on disk to the covariances file. Should be bgzipped and tabix-indexed. Note that the file need
   *    only be tabix-indexed on the start position of each row, and not the range of each row (this would be difficult
   *    since neither rvtest nor RAREMETALWORKER format include an "end position" column.)
   *
   *  Once a loader object is created, call load_region() to load statistics from a specific region into memory.
   */
  SummaryStatisticsLoader(const std::vector<std::string>& score_vec, const std::vector<std::string>& cov_vec);

  /**
   * Load a region of score statistics and covariances into memory.
   * @param chromosome Chromosome.
   * @param start Integer start position of the region.
   * @param end Integer end position of the region.
   */
  void load_region(const std::string& chromosome, uint64_t start, uint64_t stop);

  // Return the covariances.
  shared_ptr<LDQueryResult> getCovResult();

  // Return the score statistics.
  shared_ptr<ScoreStatQueryResult> getScoreResult();

  /**
   * Getter to return the residual variance under the null model.
   * @return sigma2
   */
  double getSigma2() { return sigma2; }

  /**
   * Getter to return number of samples used when calculating scores/covariances.
   * @return nsamples
   */
  uint64_t getNumSamples() { return nsamples; }
};

#endif //LDSERVER_SUMMARYSTATISTICSLOADER_H
