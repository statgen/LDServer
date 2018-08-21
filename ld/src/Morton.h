//
// Created by dtaliun on 6/20/18.
//

#ifndef LDSERVER_MORTON_H
#define LDSERVER_MORTON_H

#include <iostream>
#include <set>
#include <queue>

uint64_t split_bits(uint64_t value);

uint64_t combine_bits(uint64_t value);

uint64_t to_morton_code(uint64_t x, uint64_t y);

void from_morton_code(uint64_t z, uint64_t& x, uint64_t& y);

uint64_t load_bits(uint64_t bit_pattern, uint32_t bit_position, uint64_t value, uint32_t dim);

uint64_t compute_bigmin(uint64_t xd, uint64_t z_min, uint64_t z_max);

uint64_t compute_litmax(uint64_t xd, uint64_t z_min, uint64_t z_max);

uint64_t compute_litmax_bigmin(uint64_t xd, uint64_t z_min, uint64_t z_max, uint64_t& litmax, uint64_t& bigmin);

void get_cells(uint64_t region_start_bp, uint64_t region_stop_bp, std::set<uint64_t>& cells);

void get_cells(uint64_t index_bp, uint64_t region_start_bp, uint64_t region_end_bp, std::set<uint64_t>& cells);

#endif //LDSERVER_MORTON_H
