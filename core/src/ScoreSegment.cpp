#include "ScoreSegment.h"

ScoreSegment::ScoreSegment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp, genotypes_store store) : Segment(chromosome, start_bp, stop_bp, store) {
  score_results = make_shared<vector<ScoreResult>>();
}

ScoreSegment::ScoreSegment(Segment&& other) noexcept : Segment(std::move(other)) {
  score_results = make_shared<vector<ScoreResult>>();
}

void ScoreSegment::load(redisContext *redis_cache, const string& key) {
  redisReply *reply = nullptr;
  reply = (redisReply *) redisCommand(redis_cache, "GET %b", key.c_str(), key.length());
  if (reply == nullptr) {
    throw runtime_error("Error while reading a segment from Redis cache");
  }
  if (reply->type == REDIS_REPLY_ERROR) {
    throw runtime_error("Error while reading a segment from Redis cache: " + string(reply->str));
  }
  if (reply->len > 0) {
    stringbuf buffer(string(reply->str, reply->len), ios::binary | ios::in);
    basic_istream<char> is(&buffer);
    {
      cereal::BinaryInputArchive iarchive(is);
      load(iarchive);
    }
    cached = true;
    names_loaded = true;
    genotypes_loaded = false;
  } else {
    cached = false;
    names_loaded = false;
    genotypes_loaded = false;
  }
  freeReplyObject(reply);
}

void ScoreSegment::save(redisContext *redis_cache, const string& key) {
  redisReply *reply = nullptr;
  stringstream os(ios::binary | ios::out);
  {
    cereal::BinaryOutputArchive oarchive(os);
    save(oarchive);
  }
  string temp = os.str();
  size_t temp_size = os.tellp();
  reply = (redisReply *) redisCommand(redis_cache, "SET %b %b", key.c_str(), key.length(), temp.c_str(), temp_size);
  if (reply == nullptr) {
    throw runtime_error("Error while writing a segment to Redis cache");
  }
  if (reply->type == REDIS_REPLY_ERROR) {
    throw runtime_error("Error while writing a segment to Redis cache: " + string(reply->str));
  }
  cached = true;
  freeReplyObject(reply);
}

bool ScoreSegment::has_scores() const {
  return !score_results->empty();
}

void ScoreSegment::compute_scores(const arma::vec &phenotype) {
  // If the segment has no genotypes, we can't calculate anything.
  if (n_haplotypes == 0) {
    return;
  }

  // Load genotypes.
  arma::fmat genotypes(this->get_genotypes());

  // Center and calculate sigma2 from phenotype
  double pheno_mean = arma::mean(phenotype);
  arma::vec phenotype_centered = phenotype - pheno_mean;
  double sigma2 = arma::var(phenotype_centered, 1);

  // Calculate statistics for each variant
  for (uint64_t col = 0; col < genotypes.n_cols; col++) {
    ScoreResult result;
    result.variant = this->names[col];
    result.position = this->positions[col];
    result.chrom = this->chromosome;

    arma::vec genotype_col = arma::conv_to<arma::vec>::from(genotypes.col(col));

    // Mean center/impute
    double mean = means[col];
    genotype_col -= mean;
    genotype_col.replace(arma::datum::nan, 0);

    // Score stat
    double u = arma::dot(genotype_col, phenotype_centered);

    // Calculate denominator
    double denom = 0;
    double value;
    for (uint64_t i = 0; i < genotype_col.n_elem; i++) {
      value = genotype_col[i];
      denom += value * value;
    }
    denom = denom * sigma2;

    double v = sqrt(denom);
    double t = (u / v);
    double pvalue = 2 * arma::normcdf(-fabs(t));

    result.score_stat = u / sigma2; // match RAREMETAL convention
    result.pvalue = pvalue;
    result.alt_freq = this->freqs[col];;

    score_results->emplace_back(result);
  }
}

void ScoreSegment::add_score(ScoreResult score) {
  score_results->emplace_back(score);
}

void ScoreSegment::extract(uint64_t start, uint64_t end, struct ScoreStatQueryResult& result) const {
  int i_start, i_end;

  // If this segment doesn't have scores, or doesn't overlap the region requested,
  // signal with a sentinel value of -1 that loading stopped immediately.
  // As a side effect, overlaps_region will store the start and end index of the segment to iterate over.
  if (!this->has_scores() || !this->overlaps_region(start, end, i_start, i_end)) {
    result.last_i = -1;
    return;
  }

  int64_t i = result.last_i >= 0 ? result.last_i : i_start;
  for (; i <= i_end; i++) {
    result.data.emplace_back((*score_results)[i]);
    if (result.data.size() >= result.limit) {
      result.last_i = i + 1;
      return;
    }
  }

  result.last_i = -1;
}

bool ScoreSegment::operator==(const ScoreSegment& other) const {
  for (int i = 0; i < this->get_n_variants(); ++i) {
    if (other.get_name(i) != other.get_name(i)) { return false; }
    if (other.get_position(i) != other.get_position(i)) { return false; }
  }

  if (this->get_store() != other.get_store()) {
    return false;
  }

  double* this_value;
  double* other_value;
  double diff;
  for (int i = 0; i < this->score_results->size(); i++) {
    this_value = &(*this->score_results)[i].pvalue;
    other_value = &(*other.score_results)[i].pvalue;
    diff = fabs(*this_value - *other_value);
    if (diff > 0.00001) {
      return false;
    }
  }

  return true;
}