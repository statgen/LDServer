#include "Cell.h"

Cell::Cell(uint32_t unique_key, const string& samples_name, const string &chromosome, uint64_t morton_code) :
        key(""), cached(false), raw_fmat(nullptr), segment_i(nullptr), segment_j(nullptr), chromosome(chromosome), morton_code(morton_code) {
    stringstream os(ios::binary | ios::out);
    os.write(reinterpret_cast<const char*>(&unique_key), sizeof(unique_key));
    os.write(samples_name.c_str(), samples_name.size());
    os.write(chromosome.c_str(), chromosome.size());
    os.write(reinterpret_cast<const char*>(&morton_code), sizeof(morton_code));
    os.flush();
    key = os.str();
    from_morton_code(morton_code, this->i, this->j);
}

Cell::~Cell() {

}

const char* Cell::get_key() const {
    return key.c_str();
}

uint64_t Cell::get_key_size() const {
    return key.length();
}

void Cell::load(redisContext* redis_cache) {
    redisReply *reply = nullptr;
    reply = (redisReply *) redisCommand(redis_cache, "GET %b", key.c_str(), key.length());
    if (reply == nullptr) {
        throw runtime_error("Error while reading a cell from Redis cache");
    }
    if (reply->type == REDIS_REPLY_ERROR) {
        throw runtime_error("Error while reading a cell from Redis cache: " + string(reply->str));
    }
    if (reply->len > 0) {
        float *old_raw_mat = raw_fmat.release();
        if (old_raw_mat != nullptr) {
            delete[] old_raw_mat;
        }
        raw_fmat = unique_ptr<float[]>(new float[reply->len / sizeof(float)]);
        memcpy(reinterpret_cast<void *>(raw_fmat.get()), reply->str, reply->len);
        cached = true;
    } else {
        cached = false;
    }
    freeReplyObject(reply);
}

void Cell::save(redisContext* redis_cache) {
    redisReply* reply = nullptr;
    auto n = segment_i->positions.size() * sizeof(float);
    if (i != j) {
        n *= segment_j->positions.size();
    } else {
        n *= segment_i->positions.size();
    }
    reply = (redisReply*)redisCommand(redis_cache, "SET %b %b", key.c_str(), key.length(), raw_fmat.get(), n);
    if (reply == nullptr) {
        throw runtime_error("Error while writing a cell to Redis cache");
    }
    if (reply->type == REDIS_REPLY_ERROR) {
        throw runtime_error("Error while writing a cell to Redis cache: " + string(reply->str));
    }
    cached = true;
    freeReplyObject(reply);
}

bool Cell::is_cached() const {
    return cached;
}

void Cell::compute() {
//    auto start = std::chrono::system_clock::now();
    auto n_variants_i = segment_i->names.size();
    if (n_variants_i <= 0) {
        return;
    }
    if (this->i == this->j) { // diagonal cell
        arma::sp_fmat S_i(arma::uvec(segment_i->sp_mat_rowind.data(), segment_i->sp_mat_rowind.size(), false, false),
                          arma::uvec(segment_i->sp_mat_colind.data(), segment_i->sp_mat_colind.size(), false, false),
                          arma::fvec(segment_i->sp_mat_rowind.size(), arma::fill::ones),
                          segment_i->n_haplotypes, n_variants_i);
        arma::frowvec J(segment_i->n_haplotypes, arma::fill::ones); // vector of 1's
        arma::fmat C1(J * S_i); // allele1 counts per variant
        arma::fmat C2(segment_i->n_haplotypes - C1); // allele2 counts per variant
        arma::fmat M1(C1.t() * C1); // denominator
        arma::fmat R((segment_i->n_haplotypes * S_i.t() * S_i - M1) / sqrt(M1 % (C2.t() * C2)));
        float* old_raw_mat = raw_fmat.release();
        if (old_raw_mat != nullptr) {
            delete[] old_raw_mat;
        }
        raw_fmat = unique_ptr<float[]>(new float[R.n_elem]);
        memcpy(reinterpret_cast<void*>(raw_fmat.get()), R.memptr(), R.n_elem * sizeof(float));
    } else {
        auto n_variants_j = segment_j->names.size();
        if (n_variants_j <= 0) {
            return;
        }
        arma::sp_fmat S_i(arma::uvec(segment_i->sp_mat_rowind.data(), segment_i->sp_mat_rowind.size(), false, false),
                          arma::uvec(segment_i->sp_mat_colind.data(), segment_i->sp_mat_colind.size(), false, false),
                          arma::fvec(segment_i->sp_mat_rowind.size(), arma::fill::ones),
                          segment_i->n_haplotypes, n_variants_i);
        arma::frowvec J(segment_i->n_haplotypes, arma::fill::ones); // vector of 1's
        arma::sp_fmat S_j(
                arma::uvec(segment_j->sp_mat_rowind.data(), segment_j->sp_mat_rowind.size(), false, false),
                arma::uvec(segment_j->sp_mat_colind.data(), segment_j->sp_mat_colind.size(), false, false),
                arma::fvec(segment_j->sp_mat_rowind.size(), arma::fill::ones),
                segment_j->n_haplotypes, n_variants_j);
        arma::fmat S_i_C1(J * S_i); // allele 1 counts for segment_i variant
        arma::fmat S_i_C2(segment_i->n_haplotypes - S_i_C1); // allele 2 counts for lead variant
        arma::fmat S_j_C1(J * S_j); // allele 1 counts for segment_j variants
        arma::fmat S_j_C2(segment_j->n_haplotypes - S_j_C1); // allele 2 counts for region variants
        arma::fmat M1(S_i_C1.t() * S_j_C1);
        arma::fmat R((segment_i->n_haplotypes * S_i.t() * S_j - M1) / sqrt(M1 % (S_i_C2.t() * S_j_C2)));
        float* old_raw_mat = raw_fmat.release();
        if (old_raw_mat != nullptr) {
            delete[] old_raw_mat;
        }
        raw_fmat = unique_ptr<float[]>(new float[R.n_elem]);
        memcpy(reinterpret_cast<void*>(raw_fmat.get()), R.memptr(), R.n_elem * sizeof(float));
    }
//    auto end = std::chrono::system_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    std::cout << "Cell compute elapsed time: " << elapsed.count() << " s\n";
}

void Cell::extract(std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result) {
//    auto start = std::chrono::system_clock::now();
    if (this->i == this->j) { // diagonal cell
        if (segment_i->positions.empty()) {
            result.last_i = result.last_j = -1;
            return;
        }
        int segment_i_from = 0;
        int segment_i_to = segment_i->positions.size() - 1;
        if ((region_start_bp > segment_i->start_bp) && (region_start_bp <= segment_i->stop_bp)) {
            segment_i_from = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), region_start_bp) - segment_i->positions.begin();
        }
        if (segment_i_from >= segment_i->positions.size()) {
            result.last_i = result.last_j = -1;
            return;
        }
        if ((region_stop_bp > segment_i->start_bp) && (region_stop_bp <= segment_i->stop_bp)) {
            segment_i_to = std::upper_bound(segment_i->positions.begin(), segment_i->positions.end(), region_stop_bp) - segment_i->positions.begin() - 1;
        }
        if (segment_i_to < 0) {
            result.last_i = result.last_j = -1;
            return;
        }
        int i = result.last_i >= 0 ? result.last_i : 0;
        int j = result.last_j >= 0 ? result.last_j : i + 1;
        arma::fmat R(raw_fmat.get(), segment_i->positions.size(), segment_i->positions.size(), false, true);
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto result_i = result.data.size();
        while (i < segment_i_n_variants - 1u) {
            while (j < segment_i_n_variants) {
                result.data.emplace_back(
                        segment_i->names[segment_i_from + i],
                        chromosome,
                        segment_i->positions[segment_i_from + i],
                        segment_i->names[segment_i_from + j],
                        chromosome,
                        segment_i->positions[segment_i_from + j],
                        R(segment_i_from + i, segment_i_from + j),
                        pow(R(segment_i_from + i, segment_i_from + j), 2.0)
                );
                ++result_i;
                ++j;
                if (result_i >= result.limit) {
                    if (j < segment_i_n_variants) {
                        result.last_i = i;
                        result.last_j = j;
                    } else if (++i < segment_i_n_variants - 1u) {
                        result.last_i = i;
                        result.last_j = i + 1;
                    } else {
                        result.last_i = result.last_j = -1;
                    }
                    return;
                }
            }
            ++i;
            j = i + 1;
        }
    } else {
        if ((segment_i->positions.empty()) || (segment_j->positions.empty())) {
            result.last_i = result.last_j = -1;
            return;
        }
        int segment_i_from = 0;
        int segment_i_to = segment_i->names.size() - 1;
        if ((region_start_bp > segment_i->start_bp) && (region_start_bp <= segment_i->stop_bp)) {
            segment_i_from = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), region_start_bp) - segment_i->positions.begin();
        }
        if (segment_i_from >= segment_i->positions.size()) {
            result.last_i = result.last_j = -1;
            return;
        }
        if ((region_stop_bp > segment_i->start_bp) && (region_stop_bp <= segment_i->stop_bp)) {
            segment_i_to = std::upper_bound(segment_i->positions.begin(), segment_i->positions.end(), region_stop_bp) - segment_i->positions.begin() - 1u;
        }
        if (segment_i_to < 0) {
            result.last_i = result.last_j = -1;
            return;
        }
        int segment_j_from = 0;
        int segment_j_to = segment_j->names.size() - 1;
        if ((region_start_bp > segment_j->start_bp) && (region_start_bp <= segment_j->stop_bp)) {
            segment_j_from = std::lower_bound(segment_j->positions.begin(), segment_j->positions.end(), region_start_bp) - segment_j->positions.begin();
        }
        if (segment_j_from >= segment_j->positions.size()) {
            result.last_i = result.last_j = -1;
            return;
        }
        if ((region_stop_bp > segment_j->start_bp) && (region_stop_bp <= segment_j->stop_bp)) {
            segment_j_to = std::upper_bound(segment_j->positions.begin(), segment_j->positions.end(), region_stop_bp) - segment_j->positions.begin() - 1u;
        }
        if (segment_j_to < 0) {
            result.last_i = result.last_j = -1;
            return;
        }
        int i = result.last_i >= 0 ? result.last_i : 0;
        int j = result.last_j >= 0 ? result.last_j : 0;
        arma::fmat R(raw_fmat.get(), segment_i->positions.size(), segment_j->positions.size(), false, true);
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto segment_j_n_variants = segment_j_to - segment_j_from + 1;
        auto result_i = result.data.size();
        while (i < segment_i_n_variants) {
            while (j < segment_j_n_variants) {
                result.data.emplace_back(
                        segment_i->names[segment_i_from + i],
                        chromosome,
                        segment_i->positions[segment_i_from + i],
                        segment_j->names[segment_j_from + j],
                        chromosome,
                        segment_j->positions[segment_j_from + j],
                        R(segment_i_from + i, segment_j_from + j),
                        pow(R(segment_i_from + i, segment_j_from + j), 2.0)
                );
                ++result_i;
                ++j;
                if (result_i >= result.limit) {
                    if (j < segment_j_n_variants) {
                        result.last_i = i;
                        result.last_j = j;
                    } else if (++i < segment_i_n_variants) {
                        result.last_i = i;
                        result.last_j = 0;
                    } else {
                        result.last_i = result.last_j = -1;
                    }
                    return;
                }

            }
            ++i;
            j = 0;
        }
    }
    result.last_i = result.last_j = -1;
//    auto end = std::chrono::system_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    std::cout << "Cell region extract elapsed time: " << elapsed.count() << " s\n";
}

void Cell::extract(const std::string& index_variant, std::uint64_t index_bp, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result) {
//    auto start = std::chrono::system_clock::now();
    if (this->i == this->j) { // diagonal cell
        if (segment_i->positions.empty()) {
            result.last_j = -1;
            return;
        }
        int segment_i_from = 0;
        int segment_i_to = segment_i->names.size() - 1;
        if ((region_start_bp > segment_i->start_bp) && (region_start_bp <= segment_i->stop_bp)) {
            segment_i_from = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), region_start_bp) - segment_i->positions.begin();
        }
        if (segment_i_from >= segment_i->positions.size()) {
            result.last_j = -1;
            return;
        }
        if ((region_stop_bp > segment_i->start_bp) && (region_stop_bp <= segment_i->stop_bp)) {
            segment_i_to = std::upper_bound(segment_i->positions.begin(), segment_i->positions.end(), region_stop_bp) - segment_i->positions.begin() - 1u;
        }
        if (segment_i_to < 0) {
            result.last_j = -1;
            return;
        }
        auto index_it = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), index_bp);
        if (index_it == segment_i->positions.end()) {
            result.last_j = -1;
            return;
        }
        int segment_i_index = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), index_bp) - segment_i->positions.begin();
        while (segment_i_index < segment_i->positions.size()) {
            if (segment_i->positions[segment_i_index] == index_bp) {
                if (segment_i->names[segment_i_index].compare(index_variant) == 0) {
                    break;
                }
                ++segment_i_index;
            } else {
                segment_i_index = segment_i->positions.size();
            }
        }
        if (segment_i_index >= segment_i->positions.size()) {
            result.last_j = -1;
            return;
        }
        arma::fmat R(raw_fmat.get(), segment_i->positions.size(), segment_i->positions.size(), false, true);
        int i = segment_i_index;
        int j = result.last_j >= 0 ? result.last_j : 0;
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto result_i = result.data.size();
        while (j < segment_i_n_variants) {
            result.data.emplace_back(
                    segment_i->names[i],
                    chromosome,
                    segment_i->positions[i],
                    segment_i->names[segment_i_from + j],
                    chromosome,
                    segment_i->positions[segment_i_from + j],
                    R(i, segment_i_from + j),
                    pow(R(i, segment_i_from + j), 2.0)
            );
            ++result_i;
            ++j;
            if (result_i >= result.limit) {
                if (j < segment_i_n_variants) {
                    result.last_j = j;
                } else {
                    result.last_j = -1;
                }
                return;
            }
        }
    } else {
        if ((segment_i->positions.empty()) || (segment_j->positions.empty())) {
            result.last_j = -1;
            return;
        }

        shared_ptr<Segment> segment_i;
        shared_ptr<Segment> segment_j;

        bool reversed = false;

        if ((index_bp >= this->segment_i->start_bp) && (index_bp <= this->segment_i->stop_bp)) {
            segment_i = this->segment_i;
            segment_j = this->segment_j;
        } else if ((index_bp >= this->segment_j->start_bp) && (index_bp <= this->segment_j->stop_bp)) {
            segment_i = this->segment_j;
            segment_j = this->segment_i;
            reversed = true;
        } else {
            result.last_j = -1;
            return;
        }
        int segment_j_from = 0;
        int segment_j_to = segment_j->names.size() - 1;
        if ((region_start_bp > segment_j->start_bp) && (region_start_bp <= segment_j->stop_bp)) {
            segment_j_from = std::lower_bound(segment_j->positions.begin(), segment_j->positions.end(), region_start_bp) - segment_j->positions.begin();
        }
        if (segment_j_from >= segment_j->positions.size()) {
            result.last_j = -1;
            return;
        }
        if ((region_stop_bp > segment_j->start_bp) && (region_stop_bp <= segment_j->stop_bp)) {
            segment_j_to = std::upper_bound(segment_j->positions.begin(), segment_j->positions.end(), region_stop_bp) - segment_j->positions.begin() - 1u;
        }
        if (segment_j_to < 0) {
            result.last_j = -1;
            return;
        }
        int segment_i_index = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), index_bp) - segment_i->positions.begin();
        while (segment_i_index < segment_i->positions.size()) {
            if (segment_i->positions[segment_i_index] == index_bp) {
                if (segment_i->names[segment_i_index].compare(index_variant) == 0) {
                    break;
                }
                ++segment_i_index;
            } else {
                segment_i_index = segment_i->positions.size();
            }
        }
        if (segment_i_index >= segment_i->positions.size()) {
            result.last_j = -1;
            return;
        }
        arma::fmat R(raw_fmat.get(), this->segment_i->positions.size(), this->segment_j->positions.size(), false, true);
        int i = segment_i_index;
        int j = result.last_j >= 0 ? result.last_j : 0;
        auto segment_j_n_variants = segment_j_to - segment_j_from + 1;
        auto result_i = result.data.size();
        while (j < segment_j_n_variants) {
            result.data.emplace_back(
                    segment_i->names[i],
                    chromosome,
                    segment_i->positions[i],
                    segment_j->names[segment_j_from + j],
                    chromosome,
                    segment_j->positions[segment_j_from + j],
                    reversed ? R(segment_j_from + j, i) : R(i, segment_j_from + j),
                    pow(reversed ? R(segment_j_from + j, i) : R(i, segment_j_from + j), 2.0)
            );
            ++result_i;
            ++j;
            if (result_i >= result.limit) {
                if (j < segment_j_n_variants) {
                    result.last_j = j;
                } else {
                    result.last_j = -1;
                }
                return;
            }
        }
    }
    result.last_j = -1;
//    auto end = std::chrono::system_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    std::cout << "Cell variant extract elapsed time: " << elapsed.count() << " s\n";
}


