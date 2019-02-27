#ifndef LDSERVER_CELL_H
#define LDSERVER_CELL_H

#include <stdexcept>
#include <iostream>
#include <memory>
#include <map>
#include <string>
#include <algorithm>
#include <chrono>
#include <savvy/reader.hpp>
#include <savvy/armadillo_vector.hpp>
#include <hiredis/hiredis.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include "Raw.h"
#include "Morton.h"
#include "Segment.h"
#include "Types.h"
#include <bitset>

using namespace std;

/**
 * Class to represent a "cell" of a matrix of segments, where each segment is a fixed width chunk of the genome.
 * The indexes i and j correspond to the position in that matrix. Often these indexes are linearized into a single
 * index z by using morton codes.
 * https://en.wikipedia.org/wiki/Z-order_curve
 */
class Cell {
protected:
    bool cached;
    uint64_t i;
    uint64_t j;

    /**
     * Matrix of computed correlation values (r, r^2, or covariance, etc.)
     * If segment i has N rows and M columns, and segment j has O rows and P columns, then
     * raw_fmat will be a matrix of dimensions M x P.
     */
    unique_ptr<float[]> raw_fmat;

public:
    shared_ptr<Segment> segment_i;
    shared_ptr<Segment> segment_j;

    Cell(uint64_t i, uint64_t j);
    virtual ~Cell();

    uint64_t get_i() const;
    uint64_t get_j() const;
    bool is_diagonal() const;

    /**
     * Functions to load/save to redis cache.
     *
     * The cache key is typically created from a combination of: reference panel or genotype dataset ID,
     * name of sample subset, correlation type, chromosome, and morton code.
     *
     * Only the raw_fmat (matrix of computed correlation values) is stored in the cache.
     *
     * @param redis_cache
     * @param key
     */
    void load(redisContext* redis_cache, const string& key);
    void save(redisContext* redis_cache, const string& key);

    bool is_cached() const;

    virtual void compute() = 0;

    void extract(std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result, bool diagonal = false);
    void extract(const std::string& index_variant, std::uint64_t index_bp, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result);

};

class CellR : public Cell {
public:
    using Cell::Cell;
    virtual ~CellR();
    void compute() override;
};

class CellRsquare : public Cell {
public:
    using Cell::Cell;
    virtual ~CellRsquare();
    void compute() override;
};

class CellCov : public Cell {
public:
    using Cell::Cell;
    virtual ~CellCov();
    void compute() override;
};

class CellRsquareApprox : public Cell {
public:
    using Cell::Cell;
    virtual ~CellRsquareApprox();
    void compute() override;
};

class CellFactory {
public:
    static shared_ptr<Cell> create(correlation correlation_type, uint64_t i, uint64_t j);
};
#endif
