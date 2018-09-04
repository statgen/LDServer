#include "Segment.h"

Segment::Segment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp) : key(nullptr), key_size(0u), chromosome(chromosome), start_bp(start_bp), stop_bp(stop_bp) {
    strstreambuf buffer;
    basic_ostream<char> os(&buffer);
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

