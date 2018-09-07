#ifndef LDSERVER_CELL_H
#define LDSERVER_CELL_H

#include <iostream>
#include <memory>
#include <map>
#include <string>
#include <algorithm>
#include <strstream>
#include <savvy/reader.hpp>
#include <savvy/armadillo_vector.hpp>
#include <hiredis/hiredis.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include "Raw.h"
#include "Morton.h"
#include "Segment.h"
#include "Types.h"

using namespace std;

class Cell {
private:
    char* key;
    uint64_t key_size;
    bool cached;
    arma::fmat R;

public:
    string chromosome;
    uint64_t morton_code;
    uint64_t i;
    uint64_t j;

    shared_ptr<Segment> segment_i;
    shared_ptr<Segment> segment_j;

    Cell(uint32_t unique_key, const string& samples_name, const string& chromosome, uint64_t morton_code);
    virtual ~Cell();

    const char* get_key() const;
    uint64_t get_key_size() const;

    void load(redisContext* redis_cache);
    void save(redisContext* redis_cache) const;

    bool is_cached() const;

    void compute();

    void extract(std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result);
    void extract(const std::string& index_variant, std::uint64_t index_bp, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result);

};

#endif
