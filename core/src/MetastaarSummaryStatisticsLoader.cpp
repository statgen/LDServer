#include "MetastaarSummaryStatisticsLoader.h"
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
  extract_parquet_value(s, meta, "nrows", spstoull_uint64, pq_meta.nrows, true);
  extract_parquet_value(s, meta, "ncols", spstoull_uint64, pq_meta.ncols, true);
  extract_parquet_value(s, meta, "cov_maf_cutoff", spstod, pq_meta.cov_maf_cutoff, true);
  extract_parquet_value(s, meta, "pos_start", spstoull_uint64, pq_meta.pos_start);
  extract_parquet_value(s, meta, "pos_end", spstoull_uint64, pq_meta.pos_end);
  extract_parquet_value(s, meta, "pos_mid", spstoull_uint64, pq_meta.pos_mid, true);
  extract_parquet_value(s, meta, "region_start", spstoull_uint64, pq_meta.region_start);
  extract_parquet_value(s, meta, "region_mid", spstoull_uint64, pq_meta.region_mid);
  extract_parquet_value(s, meta, "region_end", spstoull_uint64, pq_meta.region_end);

  extract_parquet_value(s, meta, "chrom", stos, pq_meta.chrom);

  return pq_meta;
}

MetastaarSummaryStatisticsLoader::MetastaarSummaryStatisticsLoader(const std::vector<std::string>& score_vec, const std::vector<std::string>& cov_vec) {
  typedef Interval<uint64_t, MetastaarParquetMetadata> MetaInterval;

  // Create map from chromosome to list of genomic intervals covered by each score statistic file.
  // Each MetaInterval contains the start/end position of the interval covered by the file, in addition to a copy of the
  // parquet metadata for each file.
  map<string, vector<MetaInterval>> score_map;
  for (auto& f : score_vec) {
    auto meta = read_parquet_metadata(f);
    MetaInterval intv(meta.region_start, meta.region_mid, meta);
    score_map[meta.chrom].emplace_back(intv);
  }

  // Construct an interval tree for each chromosome representing the score statistic file intervals.
  for (const auto& p : score_map) {
    score_tree.emplace(make_pair(p.first, p.second));
  }

  // Same as above, only now created for each covariance matrix file.
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

  uint64_t index = 0;   // Global index over all score files
  bool enteredRegion = false;

  // Imagine each "block" as the matrix of data in each score statistic file, corresponding to a range of variants with
  // some start/end position. The block variable indexes which block we are working on. The reason is that later,
  // covariance files need to know their corresponding score statistic block, but also additionally the (n+1)th block as well.
  // We keep track of a global index where 0 is the first variant, and the last index is (nvariants - 1) where nvariants
  // is the number of variants across all score statistic files that overlap the genomic range we're interested in extracting.
  vector<uint64_t> block_ends(score_overlaps.size());
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

        /* If we've made it here, we're within the requested chr:start-stop */

        // Calculate p-value from score statistic. This isn't included by default in the MetaSTAAR single variant / score stat files.
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

  // The MetaSTAAR cov matrix is sparse; it only contains covariance where the value > 0. However, this software was written
  // originally to handle raremetal/rvtest covariance, and so per region or gene we expect a covariance for every pair of
  // variants. Therefore, we need to compute for each possible combination of variants a "baseline covariance", which is
  // simply 0 - (GtU)_i * (UtG)_j, where i and j iterate over all combinations of variants. Later, as we read in values
  // from the sparse covariance matrix, we will overwrite the baseline cov computed here.
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

        // Compute the final covariance value given GtG and GtU.
        double v = (GtG - arma::dot(row_GtU, col_GtU)) / nsamples;

        // Store it to our temporary covariance store. We needed this previously to fill in all the covariance values for
        // pairs of variants where GtG = 0.
        cov_store[make_pair(row_index, col_index)] = VariantsPair(row_sstat->variant, row_sstat->chrom, row_sstat->position, col_sstat->variant, col_sstat->chrom, col_sstat->position, v);
      }
    }
  }

  // Now we have computed every possible covariance value per pair of variants in the region, we can store them
  // to the final result object.
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