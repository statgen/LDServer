#ifndef LDSERVER_RAREMETALSUMMARYSTATISTICSLOADER_H
#define LDSERVER_RAREMETALSUMMARYSTATISTICSLOADER_H

#include "SummaryStatisticsLoader.h"

class ScoreCovColumn {
protected:
  std::string name;
  uint16_t index;
public:
  ScoreCovColumn(const std::string& name, const uint16_t& index) : name(name), index(index) {}
  inline const std::string& get_name() const { return name; }
  inline const uint16_t& get_index() const { return index; }
  inline operator int() const { return index; }
  inline operator std::string() const { return name; }
};

struct ScoreColumnSpec {
  const ScoreCovColumn colChrom;
  const ScoreCovColumn colPos;
  const ScoreCovColumn colRef;
  const ScoreCovColumn colAlt;
  const ScoreCovColumn colInformativeN;
  const ScoreCovColumn colAltFreq;
  const ScoreCovColumn colInformativeAltAc;
  const ScoreCovColumn colU;
  const ScoreCovColumn colV;
  const ScoreCovColumn colEffectAllele;
  const ScoreCovColumn colPvalue;
};

struct CovColumnSpec {
  const ScoreCovColumn colChrom;
  const ScoreCovColumn colStartPos;
  const ScoreCovColumn colPos;
  const ScoreCovColumn colCov;
};

const ScoreColumnSpec SCORE_COLUMNS_RVTEST = {
  {"CHROM", 0},
  {"POS", 1},
  {"REF", 2},
  {"ALT", 3},
  {"N_INFORMATIVE", 4},
  {"AF", 5},
  {"INFORMATIVE_ALT_AC", 6},
  {"U_STAT", 12},
  {"SQRT_V_STAT", 13},
  {"effect allele", 3},
  {"PVALUE", 15}
};

const ScoreColumnSpec SCORE_COLUMNS_RAREMETAL = {
  {"CHROM", 0},
  {"POS", 1},
  {"REF", 2},
  {"ALT", 3},
  {"N_INFORMATIVE", 4},
  {"AF", 5},
  {"INFORMATIVE_ALT_AC", 7},
  {"U_STAT", 13},
  {"SQRT_V_STAT", 14},
  {"effect allele", 3},
  {"PVALUE", 16},
};

const CovColumnSpec COV_COLUMNS_RAREMETAL = {
  {"CHROM", 0},
  {"CURRENT_POS", 1},
  {"POS", 2},
  {"COV", 3},
};

const CovColumnSpec COV_COLUMNS_RVTEST = {
  {"CHROM", 0},
  {"START_POS", 1},
  {"POS", 4},
  {"COV", 5},
};

void getNthDataLine(const std::string& filepath, std::string& out, int n);
ScoreCovFormat detectScoreCovFormat(const std::string& filepath);

template <typename T>
T extract_numeric(T func(const string&), const string& value, const ScoreCovColumn& col, const string& filepath, const string& variant) {
  try {
    return func(value);
  }
  catch (...) {
    throw LDServerGenericException(
      "Invalid value detected while parsing score statistic file"
    ).set_secret(
      boost::str(boost::format("File was: %s, offending value was '%s' in column '%s' for variant '%s'") % filepath % value % col.get_name() % variant)
    );
  }
}

/**
 * Loader for RAREMETAL or rvtest summary statistic datasets.
 */
class RaremetalSummaryStatisticsLoader : public SummaryStatisticsLoader {
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
  RaremetalSummaryStatisticsLoader(const std::vector<std::string>& score_vec, const std::vector<std::string>& cov_vec);

  /**
   * Load a region of score statistics and covariances into memory.
   * @param chromosome Chromosome.
   * @param start Integer start position of the region.
   * @param end Integer end position of the region.
   */
  void load_region(const std::string& chromosome, uint64_t start, uint64_t stop) override;

  // Return the covariances.
  shared_ptr<LDQueryResult> getCovResult() override;

  // Return the score statistics.
  shared_ptr<ScoreStatQueryResult> getScoreResult() override;

  /**
   * Getter to return the residual variance under the null model.
   * @return sigma2
   */
  double getSigma2() override { return sigma2; }

  /**
   * Getter to return number of samples used when calculating scores/covariances.
   * @return nsamples
   */
  uint64_t getNumSamples() override { return nsamples; }
};

#endif
