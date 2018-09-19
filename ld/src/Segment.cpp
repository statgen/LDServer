#include "Segment.h"

Segment::Segment(uint32_t unique_key, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp) :
        key(""), cached(false), genotypes_loaded(false), variants_loaded(false), chromosome(chromosome), start_bp(start_bp), stop_bp(stop_bp) {
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
    this->variants_loaded = segment.variants_loaded;
}

Segment::~Segment() {

}

const char* Segment::get_key() const {
    return key.c_str();
}

uint64_t Segment::get_key_size() const {
    return key.length();
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
//        strstreambuf buffer(reply->str, reply->len);
        stringbuf buffer(string(reply->str, reply->len), ios::binary | ios::in);
        basic_istream<char> is(&buffer);
        {
            cereal::BinaryInputArchive iarchive(is);
            load(iarchive);
        }
        cached = true;
        variants_loaded = true;
        genotypes_loaded = false;
    } else {
        cached = false;
        variants_loaded = false;
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