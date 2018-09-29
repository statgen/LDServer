#include "Cell.h"

Cell::Cell(uint64_t i, uint64_t j) : cached(false), i(i), j(j) {

}

Cell::~Cell() {

}

uint64_t Cell::get_i() const {
    return i;
}

uint64_t Cell::get_j() const {
    return j;
}

bool Cell::is_diagonal() const {
    return i == j;
}


void Cell::load(redisContext* redis_cache, const string& key) {
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

void Cell::save(redisContext* redis_cache, const string& key) {
    redisReply* reply = nullptr;
    auto n = segment_i->get_n_variants() * sizeof(float);
    if (i != j) {
        n *= segment_j->get_n_variants();
    } else {
        n *= segment_i->get_n_variants();
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

void Cell::extract(std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result) {
//    auto start = std::chrono::system_clock::now();
    if (this->i == this->j) { // diagonal cell
        int segment_i_from = 0, segment_i_to = 0;
        if (segment_i->is_empty() || !segment_i->overlaps_region(region_start_bp, region_stop_bp, segment_i_from, segment_i_to)) {
            result.last_i = result.last_j = -1;
            return;
        }
        int i = result.last_i >= 0 ? result.last_i : 0;
        int j = result.last_j >= 0 ? result.last_j : i + 1;
        arma::fmat R(raw_fmat.get(), segment_i->get_n_variants(), segment_i->get_n_variants(), false, true);
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto result_i = result.data.size();
        while (i < segment_i_n_variants - 1u) {
            while (j < segment_i_n_variants) {
                Segment::create_pair(*segment_i, *segment_i, segment_i_from + i, segment_i_from + j, R(segment_i_from + i, segment_i_from + j), result.data);
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
        int segment_i_from = 0, segment_i_to = 0;
        int segment_j_from = 0, segment_j_to = 0;
        if (segment_i->is_empty() || segment_j->is_empty() ||
            !segment_i->overlaps_region(region_start_bp, region_stop_bp, segment_i_from, segment_i_to) ||
            !segment_j->overlaps_region(region_start_bp, region_stop_bp, segment_j_from, segment_j_to)) {
            result.last_i = result.last_j = -1;
            return;
        }
        int i = result.last_i >= 0 ? result.last_i : 0;
        int j = result.last_j >= 0 ? result.last_j : 0;
        arma::fmat R(raw_fmat.get(), segment_i->get_n_variants(), segment_j->get_n_variants(), false, true);
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto segment_j_n_variants = segment_j_to - segment_j_from + 1;
        auto result_i = result.data.size();
        while (i < segment_i_n_variants) {
            while (j < segment_j_n_variants) {
                Segment::create_pair(*segment_i, *segment_j, segment_i_from + i, segment_j_from + j, R(segment_i_from + i, segment_j_from + j), result.data);
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
        int segment_i_from = 0, segment_i_to = 0;
        int segment_i_index = 0;
        if (segment_i->is_empty() ||
            !segment_i->overlaps_region(region_start_bp, region_stop_bp, segment_i_from, segment_i_to) ||
            !segment_i->overlaps_variant(index_variant, index_bp, segment_i_index)) {
            result.last_j = -1;
            return;
        }
        arma::fmat R(raw_fmat.get(), segment_i->get_n_variants(), segment_i->get_n_variants(), false, true);
        int i = segment_i_index;
        int j = result.last_j >= 0 ? result.last_j : 0;
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto result_i = result.data.size();
        while (j < segment_i_n_variants) {
            Segment::create_pair(*segment_i, *segment_i, i, segment_i_from + j, R(i, segment_i_from + j), result.data);
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
        if ((segment_i->is_empty()) || (segment_j->is_empty())) {
            result.last_j = -1;
            return;
        }
        shared_ptr<Segment> segment_i;
        shared_ptr<Segment> segment_j;
        bool reversed = false;
        if ((index_bp >= this->segment_i->get_start_bp()) && (index_bp <= this->segment_i->get_stop_bp())) {
            segment_i = this->segment_i;
            segment_j = this->segment_j;
        } else if ((index_bp >= this->segment_j->get_start_bp()) && (index_bp <= this->segment_j->get_stop_bp())) {
            segment_i = this->segment_j;
            segment_j = this->segment_i;
            reversed = true;
        } else {
            result.last_j = -1;
            return;
        }
        int segment_j_from = 0, segment_j_to = 0;
        int segment_i_index = 0;
        if (!segment_j->overlaps_region(region_start_bp, region_stop_bp, segment_j_from, segment_j_to) ||
            !segment_i->overlaps_variant(index_variant, index_bp, segment_i_index)) {
            result.last_j = -1;
            return;
        }
        arma::fmat R(raw_fmat.get(), this->segment_i->get_n_variants(), this->segment_j->get_n_variants(), false, true);
        int i = segment_i_index;
        int j = result.last_j >= 0 ? result.last_j : 0;
        auto segment_j_n_variants = segment_j_to - segment_j_from + 1;
        auto result_i = result.data.size();
        while (j < segment_j_n_variants) {
            Segment::create_pair(*segment_i, *segment_j, i, segment_j_from + j, reversed ? R(segment_j_from + j, i) : R(i, segment_j_from + j), result.data);
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

CellR::~CellR() {

}

void CellR::compute() {
//    auto start = std::chrono::system_clock::now();
    auto n_variants_i = segment_i->get_n_variants();
    if (n_variants_i <= 0) {
        return;
    }
    if (this->i == this->j) { // diagonal cell
        arma::sp_fmat S_i = segment_i->get_genotypes();
        arma::frowvec J(segment_i->get_n_haplotypes(), arma::fill::ones); // vector of 1's
        arma::fmat C1(J * S_i); // allele1 counts per variant
        arma::fmat C2(segment_i->get_n_haplotypes() - C1); // allele2 counts per variant
        arma::fmat M1(C1.t() * C1); // denominator
        arma::fmat R((segment_i->get_n_haplotypes() * S_i.t() * S_i - M1) / sqrt(M1 % (C2.t() * C2)));
        float* old_raw_mat = raw_fmat.release();
        if (old_raw_mat != nullptr) {
            delete[] old_raw_mat;
        }
        raw_fmat = unique_ptr<float[]>(new float[R.n_elem]);
        memcpy(reinterpret_cast<void*>(raw_fmat.get()), R.memptr(), R.n_elem * sizeof(float));
    } else {
        auto n_variants_j = segment_j->get_n_variants();
        if (n_variants_j <= 0) {
            return;
        }
        arma::sp_fmat S_i = segment_i->get_genotypes();
        arma::frowvec J(segment_i->get_n_haplotypes(), arma::fill::ones); // vector of 1's
        arma::sp_fmat S_j = segment_j->get_genotypes();
        arma::fmat S_i_C1(J * S_i); // allele 1 counts for segment_i variant
        arma::fmat S_i_C2(segment_i->get_n_haplotypes() - S_i_C1); // allele 2 counts for lead variant
        arma::fmat S_j_C1(J * S_j); // allele 1 counts for segment_j variants
        arma::fmat S_j_C2(segment_j->get_n_haplotypes() - S_j_C1); // allele 2 counts for region variants
        arma::fmat M1(S_i_C1.t() * S_j_C1);
        arma::fmat R((segment_i->get_n_haplotypes() * S_i.t() * S_j - M1) / sqrt(M1 % (S_i_C2.t() * S_j_C2)));
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

CellRsquare::~CellRsquare() {

}

void CellRsquare::compute() {
//    auto start = std::chrono::system_clock::now();
    auto n_variants_i = segment_i->get_n_variants();
    if (n_variants_i <= 0) {
        return;
    }
    if (this->i == this->j) { // diagonal cell
        arma::sp_fmat S_i = segment_i->get_genotypes();
        arma::frowvec J(segment_i->get_n_haplotypes(), arma::fill::ones); // vector of 1's
        arma::fmat C1(J * S_i); // allele1 counts per variant
        arma::fmat C2(segment_i->get_n_haplotypes() - C1); // allele2 counts per variant
        arma::fmat M1(C1.t() * C1); // denominator
        arma::fmat R(arma::pow((segment_i->get_n_haplotypes() * S_i.t() * S_i - M1) / sqrt(M1 % (C2.t() * C2)), 2.0));
        float* old_raw_mat = raw_fmat.release();
        if (old_raw_mat != nullptr) {
            delete[] old_raw_mat;
        }
        raw_fmat = unique_ptr<float[]>(new float[R.n_elem]);
        memcpy(reinterpret_cast<void*>(raw_fmat.get()), R.memptr(), R.n_elem * sizeof(float));
    } else {
        auto n_variants_j = segment_j->get_n_variants();
        if (n_variants_j <= 0) {
            return;
        }
        arma::sp_fmat S_i = segment_i->get_genotypes();
        arma::frowvec J(segment_i->get_n_haplotypes(), arma::fill::ones); // vector of 1's
        arma::sp_fmat S_j = segment_j->get_genotypes();
        arma::fmat S_i_C1(J * S_i); // allele 1 counts for segment_i variant
        arma::fmat S_i_C2(segment_i->get_n_haplotypes() - S_i_C1); // allele 2 counts for lead variant
        arma::fmat S_j_C1(J * S_j); // allele 1 counts for segment_j variants
        arma::fmat S_j_C2(segment_j->get_n_haplotypes() - S_j_C1); // allele 2 counts for region variants
        arma::fmat M1(S_i_C1.t() * S_j_C1);
        arma::fmat R(arma::pow((segment_i->get_n_haplotypes() * S_i.t() * S_j - M1) / sqrt(M1 % (S_i_C2.t() * S_j_C2)), 2.0));
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

shared_ptr<Cell> CellFactory::create(correlation correlation_type, uint64_t i, uint64_t j) {
    switch (correlation_type) {
        case LD_R: return shared_ptr<Cell>(new CellR(i, j));
        case LD_RSQUARE: return shared_ptr<Cell>(new CellRsquare(i, j));
    }
    throw runtime_error("Unknown correlation type");
}