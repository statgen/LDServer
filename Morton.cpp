//
// Created by dtaliun on 6/20/18.
//

#include "Morton.h"

uint64_t split_bits(uint64_t value) {
    value &= 0xffffffff;
    value = (value ^ (value << 16)) & 0xffff0000ffff;
    value = (value ^ (value << 8)) & 0xff00ff00ff00ff;
    value = (value ^ (value << 4)) & 0xf0f0f0f0f0f0f0f;
    value = (value ^ (value << 2)) & 0x3333333333333333;
    value = (value ^ (value << 1)) & 0x5555555555555555;
    return value;
}

uint64_t combine_bits(uint64_t value) {
    value &= 0x5555555555555555;
    value = (value ^ (value >> 1)) & 0x3333333333333333;
    value = (value ^ (value >> 2)) & 0xf0f0f0f0f0f0f0f;
    value = (value ^ (value >> 4)) & 0xff00ff00ff00ff;
    value = (value ^ (value >> 8)) & 0xffff0000ffff;
    value = (value ^ (value >> 16)) & 0xffffffff;
    return value;
}

uint64_t to_morton_code(uint64_t x, uint64_t y) {
    return split_bits(x) | (split_bits(y) << 1);
}

void from_morton_code(uint64_t z, uint64_t& x, uint64_t& y) {
    x = combine_bits(z);
    y = combine_bits(z >> 1);
}

uint64_t load_bits(uint64_t bit_pattern, uint32_t bit_position, uint64_t value, uint32_t dim) { // dim = 0 for x; dim = 1 for y
    uint64_t wipe_mask = ~(split_bits(0xffffffff >> (32u - (bit_position / 2u + 1u))) << dim);
    bit_pattern = split_bits(bit_pattern) << dim;
    return  (value & wipe_mask) | bit_pattern;
}

uint64_t compute_bigmin(uint64_t xd, uint64_t z_min, uint64_t z_max) {
    uint64_t bigmin = 0u;
    uint64_t mask = 0x8000000000000000;
    uint32_t bit_position = 63u;
    do {
        uint64_t z_min_bit = z_min & mask;
        uint64_t z_max_bit = z_max & mask;
        uint64_t xd_bit = xd & mask;

        uint32_t dim = bit_position % 2u;
        uint64_t bit_mask = 0x1 << (bit_position / 2u);

        if (xd_bit == 0u && z_min_bit == 0u && z_max_bit == 0u) {
            // do not do anything
        } else if (xd_bit == 0u && z_min_bit == 0u && z_max_bit > 0u) {
            bigmin = load_bits(bit_mask, bit_position, z_min, dim);
            z_max = load_bits(bit_mask - 1u, bit_position, z_max, dim);
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            std::cout << "BAD!" << std::endl;
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit > 0u) {
            bigmin = z_min;
            return bigmin;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit == 0u) {
            return bigmin;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit > 0u) {
            z_min = load_bits(bit_mask, bit_position, z_min, dim);
        } else if (xd_bit > 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            std::cout << "BAD!" << std::endl;
        } else if (xd_bit > 0u && z_min_bit > 0u && z_max_bit > 0u) {
            // do not do anything
        }
        bit_position -= 1u;
        mask >>= 1;
    } while (mask != 0u);
}

uint64_t compute_litmax(uint64_t xd, uint64_t z_min, uint64_t z_max) {
    uint64_t litmax = 0u;
    uint64_t mask = 0x8000000000000000;
    uint32_t bit_position = 63u;
    do {
        uint64_t z_min_bit = z_min & mask;
        uint64_t z_max_bit = z_max & mask;
        uint64_t xd_bit = xd & mask;

        uint32_t dim = bit_position % 2u;
        uint64_t bit_mask = 0x1 << (bit_position / 2u);

        if (xd_bit == 0u && z_min_bit == 0u && z_max_bit == 0u) {
            // do not do anything
        } else if (xd_bit == 0u && z_min_bit == 0u && z_max_bit > 0u) {
            z_max = load_bits(bit_mask - 1u, bit_position, z_max, dim);
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            std::cout << "BAD!" << std::endl;
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit > 0u) {
            return litmax;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit == 0u) {
            litmax = z_max;
            return litmax;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit > 0u) {
            litmax = load_bits(bit_mask - 1u, bit_position, z_max, dim);
            z_min = load_bits(bit_mask, bit_position, z_min, dim);
        } else if (xd_bit > 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            std::cout << "BAD!" << std::endl;
        } else if (xd_bit > 0u && z_min_bit > 0u && z_max_bit > 0u) {
            // do not do anything
        }
        bit_position -= 1u;
        mask >>= 1;
    } while (mask != 0u);
}

uint64_t compute_litmax_bigmin(uint64_t xd, uint64_t z_min, uint64_t z_max, uint64_t& litmax, uint64_t& bigmin) {
    uint64_t mask = 0x8000000000000000;
    uint32_t bit_position = 63u;
    do {
        uint64_t z_min_bit = z_min & mask;
        uint64_t z_max_bit = z_max & mask;
        uint64_t xd_bit = xd & mask;

        uint32_t dim = bit_position % 2u;
        uint64_t bit_mask = 0x1 << (bit_position / 2u);

        if (xd_bit == 0u && z_min_bit == 0u && z_max_bit == 0u) {
            // do not do anything
        } else if (xd_bit == 0u && z_min_bit == 0u && z_max_bit > 0u) {
            bigmin = load_bits(bit_mask, bit_position, z_min, dim);
            z_max = load_bits(bit_mask - 1u, bit_position, z_max, dim);
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            std::cout << "BAD!" << std::endl;
            // todo:: exception
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit > 0u) {
            bigmin = z_min;
            break;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit == 0u) {
            litmax = z_max;
            break;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit > 0u) {
            litmax = load_bits(bit_mask - 1u, bit_position, z_max, dim);
            z_min = load_bits(bit_mask, bit_position, z_min, dim);
        } else if (xd_bit > 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            std::cout << "BAD!" << std::endl;
            // todo:: exception
        } else if (xd_bit > 0u && z_min_bit > 0u && z_max_bit > 0u) {
            // do not do anything
        }
        bit_position -= 1u;
        mask >>= 1;
    } while (mask != 0u);
}

void get_cells(uint64_t start_bp, uint64_t end_bp, std::set<uint64_t>& range) {
    uint64_t i = start_bp / 100u;
    uint64_t j = end_bp / 100u;
    uint64_t z_min = to_morton_code(i, i);
    uint64_t z_max = to_morton_code(j, j);
    uint64_t xd, r1, r2, bigmin, litmax;

    std::queue<std::tuple<uint64_t, uint64_t>> z_ranges;

    z_ranges.push(std::make_tuple(z_min, z_max));
    while (!z_ranges.empty()) {
        z_min = std::get<0>(z_ranges.front());
        z_max = std::get<1>(z_ranges.front());
        z_ranges.pop();
        xd = z_min + (z_max - z_min) / 2;
        from_morton_code(xd, r1, r2);
        if (r1 >= i && r1 <= j && r2 >= i && r2 <= j) {
            if (r1 <= r2) { // only upper triangle of the matrix is needed
                range.insert(xd);
            }
            if (xd > z_min) {
                z_ranges.push(std::make_pair(z_min, xd - 1));
            }
            if (xd < z_max) {
                z_ranges.push(std::make_pair(xd + 1, z_max));
            }
        } else {
            compute_litmax_bigmin(xd, z_min, z_max, litmax, bigmin);
            z_ranges.push(std::make_pair(z_min, litmax));
            z_ranges.push(std::make_pair(bigmin, z_max));
        }
    }
}

void get_cells(uint64_t index_bp, uint64_t region_start_bp, uint64_t region_stop_bp, std::set<uint64_t>& cells) {
    uint64_t index_i = index_bp / 100u, region_i = region_start_bp / 100u, region_j = region_stop_bp / 100u;
    uint64_t z_min = 0, z_max = 0;
    if (index_i <= region_i) { // only upper triangle of the matrix must be used
        z_min = to_morton_code(region_i, index_i);
        z_max = to_morton_code(region_j, index_i);
    } else if ((index_i > region_i) && (index_i <region_j)) {
        z_min = to_morton_code(index_i, region_i);
        z_max = to_morton_code(region_j, index_i);
    } else {
        z_min = to_morton_code(index_i, region_i);
        z_max = to_morton_code(index_i, region_j);
    }
    do {
        cells.insert(z_min);
        if (region_i < index_i) {
            z_min = to_morton_code(index_i, ++region_i);
        } else {
            z_min = to_morton_code(++region_i, index_i);
        }
    } while (z_min <= z_max);
}