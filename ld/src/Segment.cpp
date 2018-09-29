#include "Segment.h"

Segment::Segment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp) :
        cached(false), names_loaded(false), genotypes_loaded(false), chromosome(chromosome), start_bp(start_bp), stop_bp(stop_bp), n_haplotypes(0u) {

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
    this->sp_mat_rowind = move(segment.sp_mat_rowind);
    this->sp_mat_colind = move(segment.sp_mat_colind);
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

uint32_t Segment::get_n_variants() const {
    return names.size();
}

const string& Segment::get_name(int i) const {
    return names[i];
}

uint64_t Segment::get_position(int i) const {
    return positions[i];
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

void Segment::add(savvy::site_info& anno, savvy::armadillo::sparse_vector<float>& alleles) {
    n_haplotypes = alleles.n_rows;
    if (alleles.n_nonzero > 0) {
        std::stringstream ss("");
        ss << anno.chromosome() << ":" << anno.position() << "_" << anno.ref() << "/" << anno.alt();
        names.emplace_back(ss.str());
        positions.push_back(anno.position());
        sp_mat_colind.push_back(sp_mat_rowind.size());
        for (auto it = alleles.begin(); it != alleles.end(); ++it) {
            sp_mat_rowind.push_back(it.row());
        }
    }
}

void Segment::add_name(savvy::site_info& anno, savvy::armadillo::sparse_vector<float>& alleles) {
    if (alleles.n_nonzero > 0) {
        std::stringstream ss("");
        ss << anno.chromosome() << ":" << anno.position() << "_" << anno.ref() << "/" << anno.alt();
        names.emplace_back(ss.str());
        positions.push_back(anno.position());
    }
}

void Segment::add_genotypes(savvy::armadillo::sparse_vector<float>& alleles) {
    n_haplotypes = alleles.n_rows;
    if (alleles.n_nonzero > 0) {
        sp_mat_colind.push_back(sp_mat_rowind.size());
        for (auto it = alleles.begin(); it != alleles.end(); ++it) {
            sp_mat_rowind.push_back(it.row());
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
    sp_mat_colind.push_back(sp_mat_rowind.size());
    genotypes_loaded = true;
//    cout << "Freeze genotypes: " << names.size() << " " << positions.size() << endl;
}

bool Segment::has_names() const {
    return names_loaded;
}

bool Segment::has_genotypes() const {
    return genotypes_loaded;
}

arma::sp_fmat Segment::get_genotypes() {
    return arma::sp_fmat(arma::uvec(sp_mat_rowind.data(), sp_mat_rowind.size(), false, false),
                          arma::uvec(sp_mat_colind.data(), sp_mat_colind.size(), false, false),
                          arma::fvec(sp_mat_rowind.size(), arma::fill::ones),
                          n_haplotypes, names.size()); // hope that copy elision will take care about not using copy/move constructors
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