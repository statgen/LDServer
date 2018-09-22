#include "Segment.h"

Segment::Segment(uint32_t unique_key, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp) :
        key(""), cached(false), genotypes_loaded(false), names_loaded(false), chromosome(chromosome), start_bp(start_bp), stop_bp(stop_bp) {
    stringstream os(ios::binary | ios::out);
    os.write(reinterpret_cast<const char*>(&unique_key), sizeof(unique_key));
    os.write(samples_name.c_str(), samples_name.size());
    os.write(chromosome.c_str(), chromosome.size());
    os.write(reinterpret_cast<const char*>(&start_bp), sizeof(start_bp));
    os.write(reinterpret_cast<const char*>(&stop_bp), sizeof(stop_bp));
    os.flush();
    key = os.str();
}

Segment::Segment(Segment&& segment) {
    this->key = move(segment.key);
    this->cached = segment.cached;
    this->chromosome = move(segment.chromosome);
    this->start_bp = segment.start_bp;
    this->stop_bp = segment.stop_bp;
    this->names = move(segment.names);
    this->positions = move(segment.positions);
    this->sp_mat_rowind = move(segment.sp_mat_rowind);
    this->sp_mat_colind = move(segment.sp_mat_colind);
    this->genotypes_loaded = segment.genotypes_loaded;
    this->names_loaded = segment.names_loaded;
}

Segment::~Segment() {

}

const char* Segment::get_key() const {
    return key.c_str();
}

uint64_t Segment::get_key_size() const {
    return key.length();
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

void Segment::load(redisContext *redis_cache) {
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

void Segment::save(redisContext *redis_cache) {
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