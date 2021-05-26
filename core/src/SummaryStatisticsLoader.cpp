#include "SummaryStatisticsLoader.h"
using namespace std;

/*
 * MetaSTAAR summary statistics loading implementation
 */

template <typename T>
void extract_parquet_value(const string& file, const arrow::KeyValueMetadata& meta, const string& key, T f(const string&), T& out, bool ignore_fail=false) {
  arrow::Result<string> res = meta.Get(key);
  if (res.ok()) {
    out = f(res.ValueOrDie());
  }
  else {
    if (!ignore_fail) {
      throw LDServerGenericException(
        boost::str(boost::format("Could not extract parquet metadata for key '%s' from file, see server log for detailed exception") % key)
      ).set_secret(
        boost::str(boost::format("Failed extracting metadata key '%s' from file '%s'") % key % file)
      );
    }
  }
}

MetastaarParquetMetadata read_parquet_metadata(const string& s) {
  MetastaarParquetMetadata pq_meta;
  auto reader = parquet::ParquetFileReader::OpenFile(s);
  const arrow::KeyValueMetadata& meta = *reader->metadata()->key_value_metadata();

  pq_meta.filepath = s;
  extract_parquet_value(s, meta, "nrows", spstoull, pq_meta.nrows, true);
  extract_parquet_value(s, meta, "ncols", spstoull, pq_meta.ncols, true);
  extract_parquet_value(s, meta, "cov_maf_cutoff", spstod, pq_meta.cov_maf_cutoff, true);
  extract_parquet_value(s, meta, "pos_start", spstoull, pq_meta.pos_start);
  extract_parquet_value(s, meta, "pos_end", spstoull, pq_meta.pos_end);
  extract_parquet_value(s, meta, "pos_mid", spstoull, pq_meta.pos_mid, true);
  extract_parquet_value(s, meta, "region_start", spstoull, pq_meta.region_start);
  extract_parquet_value(s, meta, "region_mid", spstoull, pq_meta.region_mid);
  extract_parquet_value(s, meta, "region_end", spstoull, pq_meta.region_end);

  extract_parquet_value(s, meta, "chrom", stos, pq_meta.chrom);

  return pq_meta;
}

MetastaarSummaryStatisticsLoader::MetastaarSummaryStatisticsLoader(const std::vector<std::string>& score_vec, const std::vector<std::string>& cov_vec) {
  typedef Interval<uint64_t, MetastaarParquetMetadata> MetaInterval;

  map<string, vector<MetaInterval>> score_map;
  for (auto& f : score_vec) {
    auto meta = read_parquet_metadata(f);
    MetaInterval intv(meta.region_start, meta.region_mid, meta);
    score_map[meta.chrom].emplace_back(intv);
  }

  for (const auto& p : score_map) {
    score_tree.emplace(make_pair(p.first, p.second));
  }

  map<string, vector<MetaInterval>> cov_map;
  for (auto& f : cov_vec) {
    auto meta = read_parquet_metadata(f);
    MetaInterval intv(meta.region_start, meta.region_mid, meta);
    cov_map[meta.chrom].emplace_back(intv);
  }

  for (const auto& p : cov_map) {
    cov_tree.emplace(make_pair(p.first, p.second));
  }

  this->score_result = make_shared<ScoreStatQueryResult>(INIT_QUERY_LIMIT);
  this->cov_result = make_shared<LDQueryResult>(INIT_QUERY_LIMIT);
}

void MetastaarSummaryStatisticsLoader::load_region(const std::string& chromosome, uint64_t start, uint64_t stop) {
  // Pull out interval tree for particular chromosome.
  const MetastaarFileIntervalTree* chrom_score_tree;
  const MetastaarFileIntervalTree* chrom_cov_tree;
  try {
    chrom_score_tree = &score_tree.at(chromosome);
    chrom_cov_tree = &cov_tree.at(chromosome);
  }
  catch (std::out_of_range& e) {
    throw LDServerGenericException(
      boost::str(boost::format("Chromosome %s not present in score stat files") % chromosome)
    );
  }

  auto sorter = [](auto &i1, auto& i2) { return i1.start < i2.start; };

  // Find metastaar cov files overlapping the given range.
  auto cov_overlaps = chrom_cov_tree->findOverlapping(start, stop);
  if (cov_overlaps.empty()) {
    throw LDServerGenericException(boost::str(boost::format("Region %s:%s-%s did not overlap any MetaSTAAR cov file") % chromosome % start % stop));
  }

  // Retrieve list of score stat files needed as well.
  auto score_overlaps = chrom_score_tree->findOverlapping(start, stop);
  if (score_overlaps.empty()) {
    throw LDServerGenericException(boost::str(boost::format("Region %s:%s-%s did not overlap any MetaSTAAR summary stat (score) file") % chromosome % start % stop));
  }
  std::sort(score_overlaps.begin(), score_overlaps.end(), sorter);

  // Because of how the segments are created, the current cov segment needs the current score segment, AND the next score
  // segment as well.
  uint64_t off_the_end = score_overlaps.back().value.region_mid + 1;
  auto one_more = chrom_score_tree->findOverlapping(off_the_end, off_the_end);
  if (!one_more.empty()) {
    score_overlaps.push_back(one_more[0]);
  }

  // Sort overlaps just to be safe, this should be a very tiny amount of data to sort anyway.
  // We already sorted score overlaps, and then added 1 onto the end, so those should be in the proper order.
  std::sort(cov_overlaps.begin(), cov_overlaps.end(), sorter);

  // Figure out total number of variants we are going to have to load from summary stats files.
  uint64_t total_rows = 0;
  for (auto& s : score_overlaps) {
    total_rows += s.value.nrows;
  }
  uint64_t min_index = 0;
  uint64_t max_index = total_rows;

  // Allocate storage.
  vector<ScoreResult*> score_stats(total_rows);
  vector<arma::fvec> GtUm(total_rows);
  vector<uint64_t> block_ends(score_overlaps.size());

  // Global index over all score files
  uint64_t index = 0;
  bool enteredRegion = false;
  for (uint64_t block = 0; block < score_overlaps.size(); block++) {
    auto score_int = score_overlaps[block];

    // For each file found, extract data.
    // Note: we need only load data where MAF > 0 and MAF > cov_maf_cutoff. Only variants matching those two criteria
    // are stored in the MetaSTAAR covariance files.
    std::shared_ptr<arrow::io::ReadableFile> score_infile;
    PARQUET_ASSIGN_OR_THROW(score_infile, arrow::io::ReadableFile::Open(score_int.value.filepath));
    parquet::StreamReader score_reader{parquet::ParquetFileReader::Open(score_infile)};
    int ncols = score_reader.num_columns();
    int64_t nrows = score_reader.num_rows();
    int ncovariates = ncols - 10;

    // Get the covariance file corresponding to this same range.
    const auto cov_overlaps = chrom_cov_tree->findOverlapping(score_int.start, score_int.stop);
    if (cov_overlaps.size() > 1) {
      throw LDServerGenericException("Multiple MetaSTAAR covariance files overlapped a region covered by one score statistic file, should be one-to-one mapping")
        .set_secret(boost::str(boost::format("Score stat file was '%s' and region %s:%s-%s") % score_int.value.filepath % chromosome % score_int.start % score_int.stop));
    }

    // MetaSTAAR cov files only store variants with MAF > 0 and MAF < maf_cutoff.
    // The summary stat / score file is a superset of those variants.
    const double& maf_cutoff = cov_overlaps[0].value.cov_maf_cutoff;

    // Temporary storage while reading each row of parquet file. We unfortunately absolutely need to pull each
    // value while reading or it will cause a parquet reader exception.
    string chrom;
    uint32_t pos;
    string ref;
    string alt;
    uint32_t alt_AC;
    uint32_t MAC;
    double maf;
    uint32_t N;
    double U;
    double V;
    double zstat;
    double pvalue;

    double tmp;
    arma::fvec GtU(ncovariates);

    while (!score_reader.eof()) {
      // Note: you must read all values, not reading all values and sending EndRow will result in an exception
      score_reader >> chrom >> pos >> ref >> alt >> alt_AC >> MAC >> maf >> N >> U >> V;
      for (int i = 0; i < ncovariates; i++) {
        score_reader >> tmp;
        GtU[i] = tmp;
      }
      score_reader >> parquet::EndRow;

      if ((maf < 0) || (maf >= maf_cutoff)) {
        // Variants with MAF < 0 or >= cutoff are not stored in the GtG file (the "cov" file.)
        continue;
      }

      if (pos >= start) {
        if (!enteredRegion) {
          min_index = index;
          enteredRegion = true;
        }

        if (pos > stop) {
          max_index = index - 1;
          break;
        }

        // Only store required data structures / statistics if we're within the range requested.
        zstat = U / sqrt(V);
        pvalue = 2 * arma::normcdf(-fabs(zstat));

        // TODO: figure out why emplace back doesn't work with just the arguments, might give slight perf increase (?)
        ScoreResult sresult = {VariantMeta(chrom, ref, alt, pos).as_epacts(), U, pvalue, alt_AC / (2.0 * N), pos, chrom};
        score_result->data.emplace_back(sresult);
        score_stats[index] = &(score_result->data.back());
        GtUm[index] = GtU;
      }
      index++;
    }

    block_ends[block] = index - 1;
    nsamples = N;
  }

//  unordered_map<pair<string, string>, double> cov_store;
//  for (uint64_t i = 0; i < score_result->data.size(); i++) {
//    for (uint64_t j = i; j < score_result->data.size(); j++) {
//      cov_store.emplace(
//        make_pair(score_result->data[i].variant, score_result->data[j].variant),
//        0.0
//      );
//    }
//  }

  map<pair<uint64_t, uint64_t>, VariantsPair> cov_store;
  for (uint64_t i = min_index; i <= max_index; i++) {
    for (uint64_t j = i; j <= max_index; j++) {
      arma::fvec& row_GtU = GtUm[i];
      arma::fvec& col_GtU = GtUm[j];

      ScoreResult* row_sstat = score_stats[i];
      ScoreResult* col_sstat = score_stats[j];

      string& row_variant = row_sstat->variant;
      string& col_variant = col_sstat->variant;

      double v = (0.0 - arma::dot(row_GtU, col_GtU)) / nsamples;
      cov_store.emplace(make_pair(i, j), VariantsPair(row_sstat->variant, row_sstat->chrom, row_sstat->position, col_sstat->variant, col_sstat->chrom, col_sstat->position, v));
    }
  }

  for (int block = 0; block < cov_overlaps.size(); block++) {
    auto cov_int = cov_overlaps[block];

    // Now we can read the GtG file and reconstruct the complete covariance matrix.
    std::shared_ptr<arrow::io::ReadableFile> cov_infile;
    PARQUET_ASSIGN_OR_THROW(cov_infile, arrow::io::ReadableFile::Open(cov_int.value.filepath));
    parquet::StreamReader cov_reader{parquet::ParquetFileReader::Open(cov_infile)};

    // If column goes past this value, we'll need to read into the next block in order to find information on that variant.
    uint64_t& nrows = cov_int.value.nrows;

    uint32_t row;
    uint32_t col;
    double GtG;
    while (!cov_reader.eof()) {
      cov_reader >> row >> col >> GtG >> parquet::EndRow;
      uint64_t row_index = (block > 0) * (block_ends[block - 1] + 1) + row;
      uint64_t col_index = (block > 0) * (block_ends[block - 1] + 1) + col;

      bool row_in_range = (row_index >= min_index) && (row_index <= max_index);
      bool col_in_range = (col_index >= min_index) && (col_index <= max_index);

      if (row_in_range && col_in_range) {
        arma::fvec& row_GtU = GtUm[row_index];
        arma::fvec& col_GtU = GtUm[col_index];

        ScoreResult* row_sstat = score_stats[row_index];
        ScoreResult* col_sstat = score_stats[col_index];

        string& row_variant = row_sstat->variant;
        string& col_variant = col_sstat->variant;

        // TODO: inefficient, would probably be faster calculating entire matrix first - revisit later
        // Can possibly construct entire final GtU * (GtU).T matrix from ptr_aux_mem constructor in arma
        // and map position in cov -> position in GtU matrix (they won't exactly match since we only load range of positions)
        double v = (GtG - arma::dot(row_GtU, col_GtU)) / nsamples;
        cov_store[make_pair(row_index, col_index)] = VariantsPair(row_sstat->variant, row_sstat->chrom, row_sstat->position, col_sstat->variant, col_sstat->chrom, col_sstat->position, v);
      }
    }
  }

  for (uint64_t i = min_index; i <= max_index; i++) {
    for (uint64_t j = i; j <= max_index; j++) {
      ScoreResult* row_sstat = score_stats[i];
      ScoreResult* col_sstat = score_stats[j];

      cov_result->data.emplace_back(cov_store.at(make_pair(i, j)));
    }
  }

  cov_result->sort_by_variant();
}

// Return the covariances.
shared_ptr<LDQueryResult> MetastaarSummaryStatisticsLoader::getCovResult() {
  return cov_result;
}

// Return the score statistics.
shared_ptr<ScoreStatQueryResult> MetastaarSummaryStatisticsLoader::getScoreResult() {
  return score_result;
}

/**
 * Getter to return the residual variance under the null model.
 * @return sigma2
 */
double MetastaarSummaryStatisticsLoader::getSigma2() {

}

/**
 * Getter to return number of samples used when calculating scores/covariances.
 * @return nsamples
 */
uint64_t MetastaarSummaryStatisticsLoader::getNumSamples() {
  return nsamples;
}

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

void RaremetalSummaryStatisticsLoader::parseHeader(const std::string& filepath) {
  Tabix tbfile(const_cast<string&>(filepath));
  string header;
  tbfile.getHeader(header);

  /**
   * Parse program name from score file.
   */
  auto regex_prog_name = regex("##ProgramName=(\\w+)");
  auto regex_detect_rvtest = regex(".*(N_INFORMATIVE\\tAF\\tINFORMATIVE_ALT_AC).*");
  auto regex_detect_raremetal = regex(".*(N_INFORMATIVE\tFOUNDER_AF\tALL_AF\tINFORMATIVE_ALT_AC).*");

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

      cov_result->data.emplace_back(row_variant, row_chrom, row_pos, variant, row_chrom, pos, v);
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