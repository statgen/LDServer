#include "RaremetalSummaryStatisticsLoader.h"
using namespace std;

/*
 * RAREMETAL summary statistics loading implementation
 */

RaremetalSummaryStatisticsLoader::RaremetalSummaryStatisticsLoader(const std::vector<std::string>& score_vec, const std::vector<std::string>& cov_vec) {
  // Create a mapping from chromosome -> score statistic file containing that chromosome.
  for (auto& f : score_vec) {
    Tabix tb(const_cast<string&>(f));
    for (auto& chrom : tb.chroms) {
      score_map[chrom] = f;
    }
  }

  // Same as above, but for covariance files instead.
  for (auto& f : cov_vec) {
    Tabix tb(const_cast<string&>(f));
    for (auto& chrom : tb.chroms) {
      cov_map[chrom] = f;
    }
  }

  this->score_result = make_shared<ScoreStatQueryResult>(INIT_QUERY_LIMIT);
  this->cov_result = make_shared<LDQueryResult>(INIT_QUERY_LIMIT);
  this->parseHeader(score_vec[0]); // Assume all score stat files have same header
}

uint64_t RaremetalSummaryStatisticsLoader::getNumberOfVariantsFromCovFile(const string& filepath, const string& region) {
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

double RaremetalSummaryStatisticsLoader::getAltFreqForPosition(uint64_t& pos) {
  try {
    return alt_freq.at(pos);
  }
  catch (std::out_of_range& e) {
    throw LDServerGenericException("Position " + to_string(pos) + " did not have alt allele frequency when loading scores/covariance matrix");
  }
}

string RaremetalSummaryStatisticsLoader::getVariantForPosition(uint64_t& pos) {
  try {
    return pos_variant.at(pos);
  }
  catch (std::out_of_range& e) {
    throw LDServerGenericException("Position " + to_string(pos) + " does not have a known variant when loading scores/covariance matrix");
  }
}

void getNthDataLine(const std::string& filepath, std::string& out, int n) {
  unique_ptr<istream> file;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
  ifstream fs(filepath, ios_base::in | ios_base::binary);

  string line;
  inbuf.push(boost::iostreams::gzip_decompressor());
  inbuf.push(fs);
  file = make_unique<istream>(&inbuf);

  string comment = "#";
  int i = 0;
  while (getline(*file, line)) {
    if (!std::equal(comment.begin(), comment.end(), line.begin())) {
      if (i == n) {
        out = line;
        break;
      }

      i += 1;
    }
  }
}

ScoreCovFormat detectScoreCovFormat(const std::string& filepath) {
  Tabix tbfile(const_cast<string&>(filepath));
  string header;
  tbfile.getHeader(header);

  /**
   * Parse program name from score file.
   */
  auto regex_prog_name = regex("##ProgramName=(\\w+)");
  auto regex_detect_rvtest = regex(".*(N_INFORMATIVE\\tAF\\tINFORMATIVE_ALT_AC).*");
  auto regex_detect_raremetal = regex(".*(N_INFORMATIVE\tFOUNDER_AF\tALL_AF\tINFORMATIVE_ALT_AC).*");

  ScoreCovFormat detected_format;
  string format;
  smatch m1;
  if (regex_search(header, m1, regex_prog_name) && m1.size() > 1) {
    format = m1.str(1);
    if (format == "Rvtests") {
      detected_format = ScoreCovFormat::RVTEST;
    }
    else if (format == "RareMetalWorker") {
      detected_format = ScoreCovFormat::RAREMETAL;
    }
    else {
      throw LDServerGenericException(
        boost::str(boost::format("Invalid program name (%s) found in header of score statistic file") % format)
      ).set_secret(
        "Score statistic file: " + filepath
      );
    }
  }
  else {
    // Apparently there was no header containing the program name. But we can still figure it out from what column
    // names are present in the file.
    smatch m2;
    if (regex_search(header, m2, regex_detect_raremetal) && m2.size() > 1) {
      detected_format = ScoreCovFormat::RAREMETAL;
    }
    else {
      // We have to read a line from the file to get the header row, because rvtest doesn't include it in the tabix header.
      string line;
      getNthDataLine(filepath, line, 0);
      smatch m3;
      if (regex_search(line, m3, regex_detect_rvtest) && m3.size() > 1) {
        detected_format = ScoreCovFormat::RVTEST;
      }
      else {
        throw LDServerGenericException("Could not determine whether file is rvtest or raremetal format")
          .set_secret("Score statistic file: " + filepath);
      }
    }
  }

  return detected_format;
}

void RaremetalSummaryStatisticsLoader::parseHeader(const std::string& filepath) {
  Tabix tbfile(const_cast<string&>(filepath));
  string header;
  tbfile.getHeader(header);

  detected_format = detectScoreCovFormat(filepath);

  /**
   * Parse sigma2 from score file header if it is available.
   */
  if (detected_format == ScoreCovFormat::RVTEST) {
    auto regex_sigma = regex("## - Sigma2\t([0-9\\.]+)");
    smatch m4;
    if (regex_search(header, m4, regex_sigma) && !m4.empty()) {
      sigma2 = stod(m4.str(1));
    }
  }
  else {
    auto regex_sigma = regex("##Sigma_e2_Hat\t(.+)");
    smatch m5;
    if (regex_search(header, m5, regex_sigma) && !m5.empty()) {
      sigma2 = stod(m5.str(1));
    }
  }

  /**
   * Parse # of samples from score file.
   * In both rvtest and raremetalworker files, this field is exactly the same.
   */
  auto regex_samples = regex("##AnalyzedSamples=(\\d+)");
  smatch m7;
  if (regex_search(header, m7, regex_samples) && !m7.empty()) {
    nsamples = stoul(m7.str(1));
  }
  else {
    // There was no sample count in the header. Use the value in the first row of the file instead.
    // Note: with summary statistic datasets, we just return this value for informational purposes, it does not actually
    // get used in a calculation.
    string line;
    getNthDataLine(filepath, line, 1);
    vector<string> tokens;
    auto field_sep = regex("[ \t]");
    copy(sregex_token_iterator(line.begin(), line.end(), field_sep, -1), sregex_token_iterator(), back_inserter(tokens));
    nsamples = stoul(tokens[4]);
  }
}

void RaremetalSummaryStatisticsLoader::load_cov(const string& chromosome, uint64_t start, uint64_t stop) {
  cov_result->erase();

  string cov_path;
  try {
    cov_path = cov_map.at(chromosome);
  }
  catch (std::out_of_range& e) {
    throw LDServerGenericException(
      boost::str(boost::format("Chromosome %s not present in covariance matrix files") % chromosome)
    );
  }

  Tabix tbfile(const_cast<string&>(cov_path));
  string region = chromosome + ":" + to_string(start) + "-" + to_string(stop);

  bool has_chrom = find(tbfile.chroms.begin(), tbfile.chroms.end(), chromosome) != tbfile.chroms.end();
  if (!has_chrom) {
    throw LDServerGenericException("Chromosome " + chromosome + " not found within covariance matrix file");
  }

  if (alt_freq.empty()) {
    throw LDServerGenericException("No alt allele frequencies available when parsing cov matrix file - did you load the scores first?");
  }

  string line;
  vector<string> tokens;

  if ((!chromosome.empty()) && (start != 0) && (stop != 0)) {
    tbfile.setRegion(region);
  }

  const CovColumnSpec* cols;
  if (detected_format == ScoreCovFormat::RVTEST) {
    cols = &COV_COLUMNS_RVTEST;
  }
  else if (detected_format == ScoreCovFormat::RAREMETAL) {
    cols = &COV_COLUMNS_RAREMETAL;
  }

  auto separator_tab = regex("\t");
  auto separator_colon = regex(":");
  auto separator_comma = regex(",");
  while (tbfile.getNextLine(line)) {
    // Split entire line
    tokens.clear();
    copy(sregex_token_iterator(line.begin(), line.end(), separator_tab, -1), sregex_token_iterator(), back_inserter(tokens));

    string row_positions = tokens.at(cols->colPos);
    string row_cov = tokens.at(cols->colCov);
    string row_chrom = tokens.at(cols->colChrom);
    string row_startpos = tokens.at(cols->colStartPos);
    string row_chrpos = string(row_chrom).append(row_startpos);
    vector<uint64_t> positions;
    vector<string> cov_matrices;
    vector<double> cov;

    // Load/split positions
    transform(
      sregex_token_iterator(row_positions.begin(), row_positions.end(), separator_comma, -1),
      sregex_token_iterator(),
      back_inserter(positions),
      [&cols, &cov_path, &row_chrpos](const string& str) {
        return extract_numeric<int>(spstoi, str, cols->colPos, cov_path, row_chrpos);
      }
    );

    // The first element in the positions array is also the position used for this entire row
    // i.e. if positions are [3, 4, 5, 6, 7] then the covariance entires correspond to (3, 3), (3, 4), (3, 5), ...
    // always position 3 in combination with the remaining positions.
    uint64_t row_pos = positions[0];

    if ((row_pos < start) || (row_pos > stop)) {
      // This whole row stores covariances with a position we aren't interested in, so we can skip it.
      continue;
    }

    string row_variant = getVariantForPosition(row_pos);

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
      [&cols, &cov_path, &row_chrpos](const string& str) {
        return extract_numeric<double>(spstod, str, cols->colCov, cov_path, row_chrpos);
      }
    );

    // Find index of the row variant in the LD results
    auto index1 = cov_result->get_variant_index(row_variant, row_chrom, row_pos).first;

    // Load covariance data
    double row_alt_freq = getAltFreqForPosition(row_pos);
    for (uint64_t j = 0; j < cov.size(); j++) {
      uint64_t pos = positions[j];

      // Since the file is tabix-indexed by start position only, it is possible that a row will contain variants
      // beyond the chrom:start-stop that we are interested in.
      if (pos > stop) {
        break;
      }

      string variant = getVariantForPosition(pos);
      double v = cov[j];
      double j_alt_freq = getAltFreqForPosition(pos);

      /**
       * The score stats file codes variant genotypes towards the alt allele. If the alt allele frequency
       * is > 0.5, that means we're not counting towards the minor (rare) allele, and we need to flip it around.
       * We don't flip when i == j because that element represents the variance of the variant itself, which is
       * invariant to which allele we code towards (but covariance is not.)
       * We also don't flip when both the i variant and j variant need to be flipped (the ^ is XOR) because it would
       * just cancel out.
       */
      if (row_pos != pos) {
        if ((row_alt_freq > 0.5) ^ (j_alt_freq > 0.5)) {
          v = -1.0 * v;
        }
      }

      auto index2 = cov_result->get_variant_index(variant, row_chrom, pos).first;
      cov_result->add_correlation(index1, index2, v);
    }
  }
}

void RaremetalSummaryStatisticsLoader::load_scores(const string& chromosome, uint64_t start, uint64_t stop) {
  score_result->erase();
  alt_freq.clear();

  // These values are valid for the lifetime of this object, since it only exists for 1 score file and 1 cov file.
  // However, when we create a new score result for a particular region (chrom/start/stop), the erase() clears
  // these values, so they need to be restored.
  score_result->sigma2 = sigma2;
  score_result->nsamples = nsamples;

  if (start <= 0) { throw LDServerGenericException("Score statistic starting position was < 0"); }
  if (stop  <= 0) { throw LDServerGenericException("Score statistic stop position was < 0"); }

  string score_path;
  try {
    score_path = score_map.at(chromosome);
  }
  catch (std::out_of_range& e) {
    throw NoVariantsInRange(
      boost::str(boost::format("Chromosome %s not present in score statistics files") % chromosome)
    );
  }

#ifndef NDEBUG
  cout << boost::str(boost::format("Loading score statistics on chromosome %s from file %s") % chromosome % score_path) << endl;
#endif

  Tabix tbfile(const_cast<string&>(score_path));
  string region = chromosome + ":" + to_string(start) + "-" + to_string(stop);

  bool has_chrom = find(tbfile.chroms.begin(), tbfile.chroms.end(), chromosome) != tbfile.chroms.end();
  if (!has_chrom) {
    throw NoVariantsInRange("Chromosome " + chromosome + " not found within score statistic file");
  }

  string line;
  vector<string> tokens;

  if ((!chromosome.empty()) && (start != 0) && (stop != 0)) {
    tbfile.setRegion(region);
  }

  const ScoreColumnSpec* cols;
  if (detected_format == ScoreCovFormat::RVTEST) {
    cols = &SCORE_COLUMNS_RVTEST;
  }
  else if (detected_format == ScoreCovFormat::RAREMETAL) {
    cols = &SCORE_COLUMNS_RAREMETAL;
  }

  uint64_t scores_read = 0;
  while (tbfile.getNextLine(line)) {
    auto separator = regex("[ \t]");
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));

    // Extract information from line
    try {
      ScoreResult result;
      result.chrom = tokens.at(cols->colChrom);
      result.position = extract_numeric<unsigned long>(spstoul, tokens.at(cols->colPos), cols->colPos, score_path, result.variant);
      string ref = tokens.at(cols->colRef);
      string alt = tokens.at(cols->colAlt);
      result.variant = VariantMeta(result.chrom, ref, alt, result.position).as_epacts();
      result.score_stat = extract_numeric<double>(spstod, tokens.at(cols->colU), cols->colU, score_path, result.variant);
      result.pvalue = extract_numeric<double>(spstod, tokens.at(cols->colPvalue), cols->colPvalue, score_path, result.variant);

      // Get allele frequency. There is an edge case of sorts in rvtest when using related samples where the
      // BLUE estimator for singletons causes a NA result for the allele frequency. However, we don't particularly
      // care what the actual allele frequency is, only the relative frequency of the ALT vs REF allele to determine
      // which is the rare allele.
      try {
        result.alt_freq = stod(tokens.at(cols->colAltFreq));
      }
      catch (...) {
        // Allele frequency was invalid. Can we use n and alt_ac instead?
        auto n = extract_numeric<double>(spstod, tokens.at(cols->colInformativeN), cols->colInformativeN, score_path, result.variant);

        // For case/control data, the field will have `all samples : case samples : control samples`. But stod() will
        // automatically only take the first value up to the `:` and that's the value we want.
        auto alt_ac = extract_numeric<double>(spstod, tokens.at(cols->colInformativeAltAc), cols->colInformativeAltAc, score_path, result.variant);
        double af = alt_ac / (2.0 * n);
        result.alt_freq = af;
      }

      alt_freq[result.position] = result.alt_freq;
      pos_variant[result.position] = result.variant;

      // Store to ScoreStatQueryResult object
      score_result->data.emplace_back(result);
    }
    catch (LDServerGenericException& e) {
      throw;
    }
    catch(...) {
      throw LDServerGenericException(
        "Invalid value detected while parsing score statistic file"
      ).set_secret(
        boost::str(boost::format("File was: %s, offending line was:\n %s") % score_path % line)
      );
    }

    scores_read++;
    tokens.clear();
  }

  if (scores_read == 0) {
    throw NoVariantsInRange(
      boost::str(boost::format("No score statistics loaded within genomic region %s:%i-%i") % chromosome % start % stop)
    );
  }
}

void RaremetalSummaryStatisticsLoader::load_region(const std::string& chromosome, uint64_t start, uint64_t stop) {
  load_scores(chromosome, start, stop);
  load_cov(chromosome, start, stop);
}

shared_ptr<LDQueryResult> RaremetalSummaryStatisticsLoader::getCovResult() {
  return cov_result;
}

shared_ptr<ScoreStatQueryResult> RaremetalSummaryStatisticsLoader::getScoreResult() {
  return score_result;
}