#include "Segment.h"

Segment::Segment(uint32_t unique_key, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp) :
        key(nullptr), key_size(0u), cached(false), genotypes_loaded(false), variants_loaded(false), chromosome(chromosome), start_bp(start_bp), stop_bp(stop_bp) {
    strstreambuf buffer;
    basic_ostream<char> os(&buffer);
    os.write(reinterpret_cast<const char*>(&unique_key), sizeof(unique_key));
    os.write(samples_name.c_str(), samples_name.size());
    os.write(chromosome.c_str(), chromosome.size());
    os.write(reinterpret_cast<const char*>(&start_bp), sizeof(start_bp));
    os.write(reinterpret_cast<const char*>(&stop_bp), sizeof(stop_bp));
    key_size = buffer.pcount();
    key = new char[key_size];
    memcpy(key, buffer.str(), key_size);
}

Segment::Segment(Segment&& segment) {
    this->key = segment.key;
    segment.key = nullptr;
    this->key_size = segment.key_size;
    this->cached = segment.cached;
    this->chromosome = move(segment.chromosome);
    this->start_bp = segment.start_bp;
    this->stop_bp = segment.stop_bp;
    this->names = move(segment.names);
    this->positions = move(segment.positions);
    this->sp_mat_rowind = move(segment.sp_mat_rowind);
    this->sp_mat_colind = move(segment.sp_mat_colind);
}

Segment::~Segment() {
    if (key != nullptr) {
        delete[] key;
        key = nullptr;
    }
}

const char* Segment::get_key() const {
    return key;
}

uint64_t Segment::get_key_size() const {
    return key_size;
}

void Segment::load(redisContext *redis_cache) {
    redisReply *reply = nullptr;
    reply = (redisReply *) redisCommand(redis_cache, "GET %b", key, key_size);
    if (reply == nullptr) {
        throw runtime_error("Error while reading a segment from Redis cache");
    }
    if (reply->type == REDIS_REPLY_ERROR) {
        throw runtime_error("Error while reading a segment from Redis cache: " + string(reply->str));
    }
    if (reply->len > 0) {
        strstreambuf buffer(reply->str, reply->len);
        basic_istream<char> is(&buffer);
        {
            cereal::BinaryInputArchive iarchive(is);
            load(iarchive);
        }
        cached = true;
    } else {
        cached = false;
    }
    freeReplyObject(reply);
}

void Segment::save(redisContext *redis_cache) {
    redisReply *reply = nullptr;
    strstreambuf buffer;
    basic_ostream<char> os(&buffer);
    {
        cereal::BinaryOutputArchive oarchive(os);
        save(oarchive);
    }
    reply = (redisReply *) redisCommand(redis_cache, "SET %b %b", key, key_size, buffer.str(), buffer.pcount());
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