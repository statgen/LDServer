#ifndef LDSERVER_CELL_H
#define LDSERVER_CELL_H

#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <savvy/reader.hpp>
#include <savvy/armadillo_vector.hpp>
#include <hiredis/hiredis.h>
#include "Raw.h"
#include "Morton.h"
#include "Segment.h"
#include "Types.h"

using namespace std;

class Cell {
private:
    arma::fmat R;
    map<std::uint64_t, Segment>::iterator segment_i_it;
    map<std::uint64_t, Segment>::iterator segment_j_it;

public:
    string chromosome;
    uint64_t morton_code;
    uint64_t i;
    uint64_t j;

    Cell(const string& chromosome, uint64_t morton_conde);
    virtual ~Cell();

    void load(const Raw* raw, const vector<string>& samples, map<uint64_t, Segment>& segments);

    void extract(std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result);
    void extract(const std::string& index_variant, std::uint64_t index_bp, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result);
};

#endif
