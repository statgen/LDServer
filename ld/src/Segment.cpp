#include "Segment.h"

Segment::Segment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp, genotypes_store store) :
        cached(false), names_loaded(false), genotypes_loaded(false),
        chromosome(chromosome), start_bp(start_bp), stop_bp(stop_bp),
        n_haplotypes(0u), store(store) {

}

Segment::Segment(Segment&& segment) {
    this->cached = segment.cached;
    this->names_loaded = segment.names_loaded;
    this->genotypes_loaded = segment.genotypes_loaded;
    this->chromosome = move(segment.chromosome);
    this->start_bp = segment.start_bp;
    this->stop_bp = segment.stop_bp;
    this->names = move(segment.names);
    this->positions = move(segment.positions);
    this->store = segment.store;
    this->sp_mat_rowind = move(segment.sp_mat_rowind);
    this->sp_mat_colind = move(segment.sp_mat_colind);
    this->sp_mat_values = move(segment.sp_mat_values);
    this->freqs = move(segment.freqs);
    this->alleles = move(segment.alleles);
    this->alt_carriers = move(segment.alt_carriers);
    this->n_haplotypes = segment.n_haplotypes;
}

Segment::~Segment() {

}

bool Segment::is_empty() const {
    return positions.empty();
}

const string& Segment::get_chromosome() const {
    return chromosome;
}

uint64_t Segment::get_start_bp() const {
    return start_bp;
}

uint64_t Segment::get_stop_bp() const {
    return stop_bp;
}

uint64_t Segment::get_n_haplotypes() const {
    return n_haplotypes;
}

uint64_t Segment::get_n_genotypes() const {
    return n_haplotypes / 2;
}

uint32_t Segment::get_n_variants() const {
    return names.size();
}

const string& Segment::get_name(int i) const {
    return names[i];
}

uint64_t Segment::get_position(int i) const {
    return positions[i];
}

genotypes_store Segment::get_store() const {
    return store;
}

void Segment::clear() {
    clear_names();
    clear_genotypes();
}

void Segment::clear_names() {
    names.clear();
    positions.clear();
    names_loaded = false;
}

void Segment::clear_genotypes() {
    sp_mat_rowind.clear();
    sp_mat_colind.clear();
    genotypes_loaded = false;
}

void Segment::add(savvy::site_info& anno, savvy::compressed_vector<float>& alleles) {
    n_haplotypes = alleles.size();
    unsigned int n_non_zero = alleles.non_zero_size();
    if (n_non_zero > 0) {
        std::stringstream ss("");
        ss << anno.chromosome() << ":" << anno.position() << "_" << anno.ref() << "/" << anno.alt();
        names.emplace_back(ss.str());
        positions.push_back(anno.position());
        auto index_data = alleles.index_data();
        switch (store) {
            case CSC_ALL_ONES:
                sp_mat_colind.emplace_back(sp_mat_rowind.size());
                for (unsigned int i = 0u; i < n_non_zero; ++i) {
                    sp_mat_rowind.emplace_back(index_data[i]);
                }
                break;
            case CSC:
                {
                // When in CSC mode, savvy returns a vector of alleles coded 0/1/2 (additive count, or savvy:fmt:ac).
                // So the number of haplotypes is 2 * number of elements in vector, since each element represents 2 alleles.
                n_haplotypes *= 2;

                sp_mat_colind.emplace_back(sp_mat_rowind.size());
                auto value_data = alleles.value_data();
                float allele_total = 0.0;
                for (unsigned int i = 0u; i < n_non_zero; ++i) {
                    sp_mat_rowind.emplace_back(index_data[i]);
                    sp_mat_values.emplace_back(value_data[i]);
                    allele_total += value_data[i];
                }
                freqs.push_back(allele_total / (float) n_haplotypes);
                }
                break;
            case BITSET:
                this->alleles.emplace_back(n_haplotypes, false);
                alt_carriers.emplace_back();
                freqs.push_back(n_non_zero / (float)n_haplotypes);
                for (unsigned int i = 0u; i < n_non_zero; ++i) {
                    this->alleles.back()[index_data[i]] = true;
                    alt_carriers.back().emplace_back(index_data[i]);
                }
                break;
        }
    }
}

void Segment::add_name(savvy::site_info& anno, savvy::compressed_vector<float>& alleles) {
    if (alleles.non_zero_size() > 0) {
        std::stringstream ss("");
        ss << anno.chromosome() << ":" << anno.position() << "_" << anno.ref() << "/" << anno.alt();
        names.emplace_back(ss.str());
        positions.push_back(anno.position());
    }
}

void Segment::add_genotypes(savvy::compressed_vector<float>& alleles) {
    n_haplotypes = alleles.size();
    unsigned int n_non_zero = alleles.non_zero_size();
    if (n_non_zero > 0) {
        auto index_data = alleles.index_data();
        switch (store) {
            case CSC_ALL_ONES:
                sp_mat_colind.emplace_back(sp_mat_rowind.size());
                for (unsigned int i = 0u; i < n_non_zero; ++i) {
                    sp_mat_rowind.emplace_back(index_data[i]);
                }
                break;
            case CSC:
                {
                n_haplotypes *= 2;
                sp_mat_colind.emplace_back(sp_mat_rowind.size());
                auto value_data = alleles.value_data();
                for (unsigned int i = 0u; i < n_non_zero; ++i) {
                    sp_mat_rowind.emplace_back(index_data[i]);
                    sp_mat_values.emplace_back(value_data[i]);
                }
                freqs.push_back(n_non_zero / (float) n_haplotypes);
                }
                break;
            case BITSET:
                this->alleles.emplace_back(n_haplotypes, false);
                alt_carriers.emplace_back();
                freqs.push_back(n_non_zero / (float)n_haplotypes);
                for (unsigned int i = 0u; i < n_non_zero; ++i) {
                    this->alleles.back()[index_data[i]] = true;
                    alt_carriers.back().emplace_back(index_data[i]);
                }
                break;
        }
    }
}

void Segment::freeze() {
    freeze_names();
    freeze_genotypes();
}

void Segment::freeze_names() {
    names_loaded = true;
}

void Segment::freeze_genotypes() {
    if ((store == CSC_ALL_ONES) || (store == CSC)) {
        sp_mat_colind.push_back(sp_mat_rowind.size());
    }
    genotypes_loaded = true;
}

bool Segment::has_names() const {
    return names_loaded;
}

bool Segment::has_genotypes() const {
    return genotypes_loaded;
}

arma::sp_fmat Segment::get_genotypes() {
    switch (store) {
        case CSC_ALL_ONES:
            return arma::sp_fmat(arma::uvec(sp_mat_rowind.data(), sp_mat_rowind.size(), false, false),
                                 arma::uvec(sp_mat_colind.data(), sp_mat_colind.size(), false, false),
                                 arma::fvec(sp_mat_rowind.size(), arma::fill::ones),
                                 n_haplotypes, names.size()); // hope that copy elision will take care about not using copy/move constructors
        case CSC:
            return arma::sp_fmat(arma::uvec(sp_mat_rowind.data(), sp_mat_rowind.size(), false, false),
                                 arma::uvec(sp_mat_colind.data(), sp_mat_colind.size(), false, false),
                                 arma::fvec(sp_mat_values.data(), sp_mat_values.size(), false, false),
                                 n_haplotypes / 2, names.size()); // hope that copy elision will take care about not using copy/move constructors
        default:
            throw runtime_error("Not supported");
    }
}

const vector<float>& Segment::get_freqs() const {
    return freqs;
}

const vector<vector<bool>>& Segment::get_alleles() const {
    return alleles;
}

const vector<vector<unsigned int>>& Segment::get_alt_carriers() const {
    return alt_carriers;
}

void Segment::create_pair(Segment& segment1, Segment& segment2, int i, int j, double value, vector<VariantsPair>& pairs) {
    pairs.emplace_back(segment1.names[i], segment1.chromosome, segment1.positions[i], segment2.names[j], segment2.chromosome, segment2.positions[j], value);
}

bool Segment::overlaps_region(uint64_t region_start_bp, uint64_t region_stop_bp, int &segment_start_index, int &segment_stop_index) const {
    segment_start_index = 0;
    segment_stop_index = positions.size() - 1;
    if ((region_start_bp > start_bp) && (region_start_bp <= stop_bp)) {
        segment_start_index = std::lower_bound(positions.begin(), positions.end(), region_start_bp) - positions.begin();
    }
    if (segment_start_index >= positions.size()) {
        return false;
    }
    if ((region_stop_bp > start_bp) && (region_stop_bp <= stop_bp)) {
        segment_stop_index = std::upper_bound(positions.begin(), positions.end(), region_stop_bp) - positions.begin() - 1;
    }
    if (segment_stop_index < 0) {
        return false;
    }
    return true;
}

bool Segment::overlaps_variant(const string& name, uint64_t bp, int& index) const {
    index = std::lower_bound(positions.begin(), positions.end(), bp) - positions.begin();
    while (index < positions.size()) {
        if (positions[index] == bp) {
            if (names[index].compare(name) == 0) {
                return true;
            }
            ++index;
        } else {
            return false;
        }
    }
    return false;
}

void Segment::load(redisContext *redis_cache, const string& key) {
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

void Segment::save(redisContext *redis_cache, const string& key) {
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

bool Segment::is_cached() const {
    return cached;
}

ScoreSegment::ScoreSegment(Segment&& other) noexcept : Segment(std::move(other)) {
    score_results = make_shared<vector<ScoreResult>>();
}

bool ScoreSegment::has_scores() const {
    return !score_results->empty();
}

void ScoreSegment::compute_scores(const arma::vec &phenotype) {
    // Find the mean of all genotype columns.
    auto genotypes = this->get_genotypes();
    arma::fmat means(arma::mean(genotypes));

    // Mean center the genotype matrix.
    for (int i = 0; i < means.n_elem; i++) {
      genotypes.col(i) -= means[i];
    }

    // Calculate score statistics for all variants.
    arma::fmat fpheno = arma::conv_to<arma::fmat>::from(phenotype);
    auto score_stats = genotypes.t() * fpheno;

    // Calculate sigma2.
    double sigma2 = arma::var(phenotype, 1);
    double sqrt_sigma = sqrt(sigma2);

    // Calculate p-values.
    double pvalue;
    for (int i = 0; i < genotypes.n_cols; i++) {
        pvalue = 2 * arma::normcdf(-fabs(score_stats[i] / sqrt_sigma));

        ScoreResult result;
        result.score_stat = score_stats[i];
        result.sigma2 = sigma2;
        result.pvalue = pvalue;

        score_results->emplace_back(result);
    }
}