#ifndef LDSERVER_MORTON_H
#define LDSERVER_MORTON_H

#include <iostream>
#include <stdexcept>
#include <set>
#include <queue>

using namespace std;

uint64_t split_bits(uint64_t value);

uint64_t combine_bits(uint64_t value);

uint64_t to_morton_code(uint64_t x, uint64_t y);

void from_morton_code(uint64_t z, uint64_t& x, uint64_t& y);

uint64_t load_bits(uint64_t bit_pattern, uint32_t bit_position, uint64_t value, uint32_t dim);

uint64_t compute_bigmin(uint64_t xd, uint64_t z_min, uint64_t z_max);

uint64_t compute_litmax(uint64_t xd, uint64_t z_min, uint64_t z_max);

void compute_litmax_bigmin(uint64_t xd, uint64_t z_min, uint64_t z_max, uint64_t& litmax, uint64_t& bigmin);

uint64_t get_next_z(uint64_t range_start, uint64_t range_end, uint64_t z_min, uint64_t z_max, uint64_t z_init);

uint64_t get_next_z(uint64_t index, uint64_t range_start, uint64_t range_end, uint64_t z_min, uint64_t z_max, uint64_t z_init);

#endif //LDSERVER_MORTON_H
