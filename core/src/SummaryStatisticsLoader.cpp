#include "SummaryStatisticsLoader.h"
using namespace std;

SummaryStatisticsLoader::SummaryStatisticsLoader(const string& score_path, const string& cov_path) {
  this->score_path = score_path;
  this->cov_path = cov_path;
  this->score_result = make_shared<ScoreStatQueryResult>(INIT_QUERY_LIMIT);
  this->cov_result = make_shared<LDQueryResult>(INIT_QUERY_LIMIT);
}

uint64_t SummaryStatisticsLoader::getNumberOfVariantsFromCovFile(const string& filepath, const string& region) {
  Tabix tbfile(const_cast<string&>(filepath));
  tbfile.setRegion(const_cast<string&>(region));

  string line;
  vector<string> tokens;
  set<unsigned int> positions;

  while (tbfile.getNextLine(line)) {
    if (line.substr(0, 1) == "#") { continue; }
    if (line.empty()) { continue; }

    auto separator = regex("[ \t]");
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));

    auto comma = regex(",");
    string pos_str = tokens.at(4);

    transform(
      sregex_token_iterator(pos_str.begin(), pos_str.end(), comma, -1),
      sregex_token_iterator(),
      inserter(positions, positions.begin()),
      [](const auto& x) {
        return stoi(x);
      }
    );

    line.clear();
    tokens.clear();
  }

  return positions.size();
}

double SummaryStatisticsLoader::getAltFreqForPosition(uint64_t& pos) {
  try {
    return alt_freq.at(pos);
  }
  catch (std::out_of_range& e) {
    throw ScoreCovParseException("Position " + to_string(pos) + " did not have alt allele frequency when loading scores/covariance matrix");
  }
}

string SummaryStatisticsLoader::getVariantForPosition(uint64_t& pos) {
  try {
    return pos_variant.at(pos);
  }
  catch (std::out_of_range& e) {
    throw ScoreCovParseException("Position " + to_string(pos) + " does not have a known variant when loading scores/covariance matrix");
  }
}

ScoreCovFormat SummaryStatisticsLoader::detectFormat(const std::string& filepath) {
  Tabix tbfile(const_cast<string&>(filepath));
  string header;
  tbfile.getHeader(header);

  auto regex_prog_name = regex("##ProgramName=(\\w+)");
  smatch match;
  string format;
  if (regex_search(header, match, regex_prog_name) && match.size() > 1) {
    format = match.str(1);
  }
  else {
    throw runtime_error("Could not find ProgramName=XX header in file: " + filepath);
  }

  if (format == "Rvtests") {
    return ScoreCovFormat::RVTEST;
  }
  else if (format == "RareMetalWorker") {
    return ScoreCovFormat::RAREMETAL;
  }
  else {
    throw runtime_error(
      boost::str(boost::format("Invalid program name (%s) found in header of %s") % format % filepath)
    );
  }
}

void SummaryStatisticsLoader::load_cov(const string& chromosome, uint64_t start, uint64_t stop) {
  cov_result->erase();

  Tabix tbfile(const_cast<string&>(cov_path));
  ScoreCovFormat format = SummaryStatisticsLoader::detectFormat(cov_path);
  string region = chromosome + ":" + to_string(start) + "-" + to_string(stop);

  bool has_chrom = find(tbfile.chroms.begin(), tbfile.chroms.end(), chromosome) != tbfile.chroms.end();
  if (!has_chrom) {
    throw std::range_error("Chromosome " + chromosome + " not found within covariance matrix file");
  }

  string line;
  vector<string> tokens;

  if ((!chromosome.empty()) && (start != 0) && (stop != 0)) {
    tbfile.setRegion(region);
  }

  CovarianceColumns cols;
  if (format == ScoreCovFormat::RVTEST) {
    cols = COV_COLUMNS_RVTEST;
  }
  else if (format == ScoreCovFormat::RAREMETAL) {
    cols = COV_COLUMNS_RAREMETAL;
  }

  auto separator_tab = regex("\t");
  auto separator_colon = regex(":");
  auto separator_comma = regex(",");
  while (tbfile.getNextLine(line)) {
    // Split entire line
    copy(sregex_token_iterator(line.begin(), line.end(), separator_tab, -1), sregex_token_iterator(), back_inserter(tokens));

    string row_positions = tokens.at(cols.colPos);
    string row_cov = tokens.at(cols.colCov);
    string row_chrom = tokens.at(cols.colChrom);
    vector<uint64_t> positions;
    vector<string> cov_matrices;
    vector<double> cov;

    // Load/split positions
    transform(
      sregex_token_iterator(row_positions.begin(), row_positions.end(), separator_comma, -1),
      sregex_token_iterator(),
      back_inserter(positions),
      [](const string& str) { return stoi(str); }
    );

    // The first element in the positions array is also the position used for this entire row
    // i.e. if positions are [3, 4, 5, 6, 7] then the covariance entires correspond to (3, 3), (3, 4), (3, 5), ...
    // always position 3 in combination with the remaining positions.
    uint64_t row_pos = positions[0];

    // Load the covariance values on this row
    // Note: rvtest will put 3 elements within the covariance column if a binary trait was used
    // The first element is essentially G.T * G
    // The second element (G.T * C) and third element (C.T * C) are not currently used
    copy(
      sregex_token_iterator(row_cov.begin(), row_cov.end(), separator_colon, -1),
      sregex_token_iterator(),
      back_inserter(cov_matrices)
    );

    // Retrieve G.T * G
    transform(
      sregex_token_iterator(cov_matrices[0].begin(), cov_matrices[0].end(), separator_comma, -1),
      sregex_token_iterator(),
      back_inserter(cov),
      [](const string& str) { return stod(str); }
    );

    // Load covariance data
    double row_alt_freq = getAltFreqForPosition(row_pos);
    string row_variant = getVariantForPosition(row_pos);
    for (uint64_t j = 0; j < cov.size(); j++) {
      uint64_t pos = positions[j];

      if (pos > stop) {
        break;
      }

      string variant = getVariantForPosition(pos);
      double v = cov[j];
      double j_alt_freq = getAltFreqForPosition(pos);

      if (row_pos != pos) {
        if ((row_alt_freq > 0.5) ^ (j_alt_freq > 0.5)) {
          v = -1.0 * v;
        }
      }

      // VariantsPair(const string& variant1, const string& chromosome1, uint64_t position1, const string& variant2,
      // const string& chromosome2, uint64_t position2, double value):
      cov_result->data.emplace_back(row_variant, row_chrom, row_pos, variant, row_chrom, pos, v);
    }

    tokens.clear();
  }
}

void SummaryStatisticsLoader::load_scores(const string& chromosome, uint64_t start, uint64_t stop) {
  score_result->erase();
  alt_freq.clear();

  if (start <= 0) { throw std::invalid_argument("Score statistic starting position was < 0"); }
  if (stop  <= 0) { throw std::invalid_argument("Score statistic stop position was < 0"); }

  Tabix tbfile(const_cast<string&>(score_path));
  ScoreCovFormat format = SummaryStatisticsLoader::detectFormat(score_path);
  string region = chromosome + ":" + to_string(start) + "-" + to_string(stop);

  bool has_chrom = find(tbfile.chroms.begin(), tbfile.chroms.end(), chromosome) != tbfile.chroms.end();
  if (!has_chrom) {
    throw std::range_error("Chromosome " + chromosome + " not found within score statistic file");
  }

  string line;
  vector<string> tokens;

  if ((!chromosome.empty()) && (start != 0) && (stop != 0)) {
    tbfile.setRegion(region);
  }

  ScoreColumns cols;
  if (format == ScoreCovFormat::RVTEST) {
    cols = SCORE_COLUMNS_RVTEST;
  }
  else if (format == ScoreCovFormat::RAREMETAL) {
    cols = SCORE_COLUMNS_RAREMETAL;
  }

  uint64_t scores_read = 0;
  while (tbfile.getNextLine(line)) {
    auto separator = regex("[ \t]");
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));

    // Extract information from line
    ScoreResult result;
    result.chrom = tokens.at(cols.colChrom);
    result.position = stoul(tokens.at(cols.colPos));
    result.score_stat = stod(tokens.at(cols.colU));
    result.pvalue = stod(tokens.at(cols.colPvalue));
    result.alt_freq = stod(tokens.at(cols.colAltFreq));
    string ref = tokens.at(cols.colRef);
    string alt = tokens.at(cols.colAlt);
    result.variant = result.chrom + ":" + to_string(result.position) + "_" + ref + "/" + alt;
    alt_freq[result.position] = result.alt_freq;
    pos_variant[result.position] = result.variant;

    // Store to ScoreStatQueryResult object
    score_result->data.emplace_back(result);

    scores_read++;
    tokens.clear();
  }

  if (scores_read == 0) {
    throw std::range_error(
      boost::str(boost::format("No score statistics loaded within genomic region %s:%i-%i") % chromosome % start % stop)
    );
  }
}

void SummaryStatisticsLoader::load_region(const std::string& chromosome, uint64_t start, uint64_t stop) {
  load_scores(chromosome, start, stop);
  load_cov(chromosome, start, stop);
}

shared_ptr<LDQueryResult> SummaryStatisticsLoader::getCovResult() {
  return cov_result;
}

shared_ptr<ScoreStatQueryResult> SummaryStatisticsLoader::getScoreResult() {
  return score_result;
}