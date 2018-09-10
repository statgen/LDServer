#include "Cell.h"

Cell::Cell(uint32_t unique_key, const string& samples_name, const string &chromosome, uint64_t morton_code) :
        key(nullptr), key_size(0u), cached(false), segment_i(nullptr), segment_j(nullptr), chromosome(chromosome), morton_code(morton_code) {
    strstreambuf buffer;
    basic_ostream<char> os(&buffer);
    os.write(reinterpret_cast<const char*>(&unique_key), sizeof(unique_key));
    os.write(samples_name.c_str(), samples_name.size());
    os.write(chromosome.c_str(), chromosome.size());
    os.write(reinterpret_cast<const char*>(&morton_code), sizeof(morton_code));
    key_size = buffer.pcount();
    key = new char[key_size];
    memcpy(key, buffer.str(), key_size);
    from_morton_code(morton_code, this->i, this->j);
}

Cell::~Cell() {
    if (key != nullptr) {
        free(key);
        key = nullptr;
    }
}

const char* Cell::get_key() const {
    return key;
}

uint64_t Cell::get_key_size() const {
    return key_size;
}

void Cell::load(redisContext* redis_cache) {
    redisReply *reply = nullptr;
    reply = (redisReply *) redisCommand(redis_cache, "GET %b", key, key_size);
    if (reply == nullptr) {
        return; // todo: throw exception, also check reply->type == REDIS_REPLY_NIL
    }
    if (reply->len > 0) {
        strstreambuf buffer(reply->str, reply->len);
        basic_istream<char> is(&buffer);
        R.load(is, arma::arma_binary);
        cached = true;
    }
    freeReplyObject(reply);
}

void Cell::save(redisContext* redis_cache) const {
    redisReply* reply = nullptr;
    strstreambuf buffer;
    basic_ostream<char> os(&buffer);
    R.save(os, arma::arma_binary);
    reply = (redisReply*)redisCommand(redis_cache, "SET %b %b", key, key_size, buffer.str(), buffer.pcount());
    if (reply == nullptr) {
        return; // todo: throw exception
    }
    // todo: check reply->type and reply->str
    freeReplyObject(reply);
}

bool Cell::is_cached() const {
    return cached;
}

void Cell::compute() {
    auto n_variants_i = segment_i->names.size();
    arma::sp_fmat S_i(arma::uvec(segment_i->sp_mat_rowind.data(), segment_i->sp_mat_rowind.size(), false, false),
                      arma::uvec(segment_i->sp_mat_colind.data(), segment_i->sp_mat_colind.size(), false, false),
                      arma::fvec(segment_i->sp_mat_rowind.size(), arma::fill::ones),
                      segment_i->n_haplotypes, n_variants_i);
    arma::frowvec J(segment_i->n_haplotypes, arma::fill::ones); // vector of 1's
    if (this->i == this->j) { // diagonal cell
        arma::fmat C1(J * S_i); // allele1 counts per variant
        arma::fmat C2(segment_i->n_haplotypes - C1); // allele2 counts per variant
        arma::fmat M1(C1.t() * C1); // denominator
        this->R = ((segment_i->n_haplotypes * S_i.t() * S_i - M1) / sqrt(M1 % (C2.t() * C2)));
    } else {
        auto n_variants_j = segment_j->names.size();
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
        this->R = ((segment_i->n_haplotypes * S_i.t() * S_j - M1) / sqrt(M1 % (S_i_C2.t() * S_j_C2)));
    }
}

void Cell::extract(std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result) {
    if (this->i == this->j) { // diagonal cell
        int segment_i_from = 0;
        int segment_i_to = segment_i->positions.size() - 1;
        if ((region_start_bp > segment_i->start_bp) && (region_start_bp <= segment_i->stop_bp)) {
            auto start_it = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), region_start_bp);
            segment_i_from = start_it - segment_i->positions.begin();
        }
        if ((region_stop_bp > segment_i->start_bp) && (region_stop_bp <= segment_i->stop_bp)) {
            auto stop_it = std::upper_bound(segment_i->positions.begin(), segment_i->positions.end(), region_stop_bp);
            segment_i_to = stop_it - segment_i->positions.begin() - 1;
        }
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto result_i = result.data.size();
        auto result_new_size = result_i + (segment_i_n_variants * (segment_i_n_variants - 1u)) / 2u;
        if (result.last_i >= 0 && result.last_j >= 0) {
            result_new_size -= (result.last_i * (2 * (segment_i_n_variants - 1) - (result.last_i - 1))) / 2u - (result.last_j - result.last_i);
        }
        if (result_new_size == 0) {
            result.last_i = result.last_j = -1;
            return;
        } else if (result_new_size > result.limit) {
            result_new_size = result.limit;
        }
        result.data.resize(result_new_size);
        int i = result.last_i >= 0 ? result.last_i : 0;
        int j = result.last_j >= 0 ? result.last_j + 1 : i + 1;
        while (i < segment_i_n_variants - 1u) {
            while (j < segment_i_n_variants) {
                result.data[result_i].variant1 = segment_i->names[segment_i_from + i];
                result.data[result_i].chromosome1 = chromosome;
                result.data[result_i].position1 = segment_i->positions[segment_i_from + i];
                result.data[result_i].variant2 = segment_i->names[segment_i_from + j];
                result.data[result_i].chromosome2 = chromosome;
                result.data[result_i].position2 = segment_i->positions[segment_i_from + j];
                result.data[result_i].r = R(segment_i_from + i, segment_i_from + j);
                result.data[result_i].rsquare = pow(result.data[result_i].r, 2.0);
                ++result_i;
                if (result_i >= result.limit) {
                    result.last_cell = morton_code;
                    result.last_i = i;
                    result.last_j = j;
                    return;
                }
                ++j;
            }
            ++i;
            j = i + 1;
        }
    } else {
        int segment_i_from = 0;
        int segment_i_to = segment_i->names.size() - 1;
        int segment_j_from = 0;
        int segment_j_to = segment_j->names.size() - 1;
        if ((region_start_bp > segment_i->start_bp) && (region_start_bp <= segment_i->stop_bp)) {
            auto start_it = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), region_start_bp);
            segment_i_from = start_it - segment_i->positions.begin();
        }
        if ((region_stop_bp > segment_i->start_bp) && (region_stop_bp <= segment_i->stop_bp)) {
            auto stop_it = std::upper_bound(segment_i->positions.begin(), segment_i->positions.end(), region_stop_bp);
            segment_i_to = stop_it - segment_i->positions.begin() - 1u;
        }
        if ((region_start_bp > segment_j->start_bp) && (region_start_bp <= segment_j->stop_bp)) {
            auto start_it = std::lower_bound(segment_j->positions.begin(), segment_j->positions.end(), region_start_bp);
            segment_j_from = start_it - segment_j->positions.begin();
        }
        if ((region_stop_bp > segment_j->start_bp) && (region_stop_bp <= segment_j->stop_bp)) {
            auto stop_it = std::upper_bound(segment_j->positions.begin(), segment_j->positions.end(), region_stop_bp);
            segment_j_to = stop_it - segment_j->positions.begin() - 1u;
        }
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto segment_j_n_variants = segment_j_to - segment_j_from + 1;
        auto result_i = result.data.size();
        auto result_new_size = result_i + segment_i_n_variants * segment_j_n_variants;
        if ((result.last_i >= 0) && (result.last_j >= 0)) {
            result_new_size = result_new_size - result.last_i * segment_j_n_variants - (result.last_j + 1);
        }
        if (result_new_size == 0) {
            result.last_i = result.last_j = -1;
            return;
        } else if (result_new_size > result.limit) {
            result_new_size = result.limit;
        }
        result.data.resize(result_new_size);
        int i = result.last_i >= 0 ? result.last_i : 0;
        int j = result.last_j >= 0 ? result.last_j + 1 : 0;
        while (i < segment_i_n_variants) {
            while (j < segment_j_n_variants) {
                result.data[result_i].variant1 = segment_i->names[segment_i_from + i];
                result.data[result_i].chromosome1 = chromosome;
                result.data[result_i].position1 = segment_i->positions[segment_i_from + i];
                result.data[result_i].variant2 = segment_j->names[segment_j_from + j];
                result.data[result_i].chromosome2 = chromosome;
                result.data[result_i].position2 = segment_j->positions[segment_j_from + j];
                result.data[result_i].r = R(segment_i_from + i, segment_j_from + j);
                result.data[result_i].rsquare = pow(result.data[result_i].r, 2.0);
                ++result_i;
                if (result_i >= result.limit) {
                    result.last_cell = morton_code;
                    result.last_i = i;
                    result.last_j = j;
                    return;
                }
                ++j;
            }
            ++i;
            j = 0;
        }
    }
    result.last_i = result.last_j = -1;
}

void Cell::extract(const std::string& index_variant, std::uint64_t index_bp, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result) {
    if (this->i == this->j) { // diagonal cell
        int segment_i_from = 0;
        int segment_i_to = segment_i->names.size() - 1;
        if ((region_start_bp > segment_i->start_bp) && (region_start_bp <= segment_i->stop_bp)) {
            auto start_it = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), region_start_bp);
            segment_i_from = start_it - segment_i->positions.begin();
        }
        if ((region_stop_bp > segment_i->start_bp) && (region_stop_bp <= segment_i->stop_bp)) {
            auto stop_it = std::upper_bound(segment_i->positions.begin(), segment_i->positions.end(), region_stop_bp);
            segment_i_to = stop_it - segment_i->positions.begin() - 1u;
        }
        auto index_it = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), index_bp);
        if (index_it == segment_i->positions.end()) {
            // TODO: throw exception - can't find index
        }
        int segment_i_index = index_it - segment_i->positions.begin();
        while (segment_i_index < segment_i->positions.size()) {
            if (segment_i->positions[segment_i_index] == index_bp) { // TODO: must check index name as well
                break;
            }
            ++segment_i_index;
        }
        if (segment_i_index >= segment_i->positions.size()) {
            // TODO: throw exception - can't find index
        }
        auto segment_i_n_variants = segment_i_to - segment_i_from + 1;
        auto result_i = result.data.size();
        auto result_new_size = result_i + segment_i_n_variants;
        if (result.last_j >= 0) {
            result_new_size = result_new_size - (result.last_j + 1);
        }
        if (result_new_size == 0) {
            result.last_i = result.last_j = -1;
            return;
        } else if (result_new_size > result.limit) {
            result_new_size = result.limit;
        }
        result.data.resize(result_new_size);
        int i = segment_i_index;
        int j = result.last_j >= 0 ? result.last_j : 0;
        while (j < segment_i_n_variants) {
            result.data[result_i].variant1 = segment_i->names[i];
            result.data[result_i].chromosome1 = chromosome;
            result.data[result_i].position1 = segment_i->positions[i];
            result.data[result_i].variant2 = segment_i->names[segment_i_from + j];
            result.data[result_i].chromosome2 = chromosome;
            result.data[result_i].position2 = segment_i->positions[segment_i_from + j];
            result.data[result_i].r = R(i, segment_i_from + j);
            result.data[result_i].rsquare = pow(result.data[result_i].r, 2.0);
            ++result_i;
            if (result_i >= result.limit) {
                result.last_cell = morton_code;
                result.last_j = j;
                return;
            }
            ++j;
        }
    } else {
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
            // TODO: throw exception
        }
        int segment_j_from = 0;
        int segment_j_to = segment_j->names.size() - 1;
        auto start_it = std::lower_bound(segment_j->positions.begin(), segment_j->positions.end(), region_start_bp);
        segment_j_from = start_it - segment_j->positions.begin();
        auto stop_it = std::upper_bound(segment_j->positions.begin(), segment_j->positions.end(), region_stop_bp);
        segment_j_to = stop_it - segment_j->positions.begin() - 1u;
        auto index_it = std::lower_bound(segment_i->positions.begin(), segment_i->positions.end(), index_bp);
        if (index_it == segment_i->positions.end()) {
            // TODO: throw exception - can't find index
        }
        int segment_i_index = index_it - segment_i->positions.begin();
        while (segment_i_index < segment_i->positions.size()) {
            if (segment_i->positions[segment_i_index] == index_bp) { // TODO: must check index name as well
                break;
            }
            ++segment_i_index;
        }
        if (segment_i_index >= segment_i->positions.size()) {
            // TODO: throw exception - can't find index
        }
        auto segment_j_n_variants = segment_j_to - segment_j_from + 1;
        auto result_i = result.data.size();
        auto result_new_size = result_i + segment_j_n_variants;
        if (result.last_j >= 0) {
            result_new_size = result_new_size - (result.last_j + 1);
        }
        if (result_new_size == 0) {
            result.last_i = result.last_j = -1;
            return;
        } else if (result_new_size > result.limit) {
            result_new_size = result.limit;
        }
        result.data.resize(result_new_size);
        int i = segment_i_index;
        int j = result.last_j >= 0 ? result.last_j : 0;
        while (j < segment_j_n_variants) {
            result.data[result_i].variant1 = segment_i->names[i];
            result.data[result_i].chromosome1 = chromosome;
            result.data[result_i].position1 = segment_i->positions[i];
            result.data[result_i].variant2 = segment_j->names[segment_j_from + j];
            result.data[result_i].chromosome2 = chromosome;
            result.data[result_i].position2 = segment_j->positions[segment_j_from + j];
            if (reversed) {
                result.data[result_i].r = R(segment_j_from + j, i);
            } else {
                result.data[result_i].r = R(i, segment_j_from + j);
            }
            result.data[result_i].rsquare = pow(result.data[result_i].r, 2.0);
            ++result_i;
            if (result_i >= result.limit) {
                result.last_cell = morton_code;
                result.last_j = j;
                return;
            }
            ++j;
        }
    }
    result.last_i = result.last_j = -1;
}

