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
        if (xd_bit == 0u && z_min_bit == 0u && z_max_bit > 0u) {
            bigmin = load_bits(bit_mask, bit_position, z_min, dim);
            z_max = load_bits(bit_mask - 1u, bit_position, z_max, dim);
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            throw logic_error("Error while computing BIGMIN");
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit > 0u) {
            bigmin = z_min;
            return bigmin;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit == 0u) {
            return bigmin;
        } else if (xd_bit > 0u && z_min_bit == 0u && z_max_bit > 0u) {
            z_min = load_bits(bit_mask, bit_position, z_min, dim);
        } else if (xd_bit > 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            throw logic_error("Error while computing BIGMIN");
        }
        --bit_position;
        mask >>= 1;
    } while (mask != 0u);
    return bigmin;
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
        if (xd_bit == 0u && z_min_bit == 0u && z_max_bit > 0u) {
            z_max = load_bits(bit_mask - 1u, bit_position, z_max, dim);
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            throw logic_error("Error while computing LITMAX");
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
            throw logic_error("Error while computing LITMAX");
        }
        --bit_position;
        mask >>= 1;
    } while (mask != 0u);
    return litmax;
}

void compute_litmax_bigmin(uint64_t xd, uint64_t z_min, uint64_t z_max, uint64_t& litmax, uint64_t& bigmin) {
    uint64_t mask = 0x8000000000000000;
    uint32_t bit_position = 63u;
    do {
        uint64_t z_min_bit = z_min & mask;
        uint64_t z_max_bit = z_max & mask;
        uint64_t xd_bit = xd & mask;
        uint32_t dim = bit_position % 2u;
        uint64_t bit_mask = 0x1 << (bit_position / 2u);
        if (xd_bit == 0u && z_min_bit == 0u && z_max_bit > 0u) {
            bigmin = load_bits(bit_mask, bit_position, z_min, dim);
            z_max = load_bits(bit_mask - 1u, bit_position, z_max, dim);
        } else if (xd_bit == 0u && z_min_bit > 0u && z_max_bit == 0u) {
            // not possible because min <= max
            throw logic_error("Error while computing LITMAX and BIGMIN");
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
            throw logic_error("Error while computing LITMAX and BIGMIN");
        }
        --bit_position;
        mask >>= 1;
    } while (mask != 0u);
}

uint64_t get_next_z(uint64_t range_start, uint64_t range_end, uint64_t z_min, uint64_t z_max, uint64_t z_init) {
    uint64_t xd_start = 0u, xd_end = 0u;
    uint64_t xd = z_init;
    while (xd <= z_max) {
        from_morton_code(xd, xd_start, xd_end);
        if (xd_start >= range_start && xd_start <= range_end && xd_end >= range_start && xd_end <= range_end) {
            if (xd_start <= xd_end) { // only upper triangle of the matrix is needed
                return xd;
            }
            ++xd;
        } else {
            xd = compute_bigmin(xd, z_min, z_max);
        }
    }
    return xd;
}

uint64_t get_next_z(uint64_t index, uint64_t range_start, uint64_t range_end, uint64_t z_min, uint64_t z_max, uint64_t z_init) {
    uint64_t xd_start = 0u, xd_end = 0u;
    uint64_t xd = z_init;
    while (xd <= z_max) {
        from_morton_code(xd, xd_start, xd_end);
        if (index <= range_start) {
            if ((xd_start >= range_start) && (xd_start <= range_end) && (index == xd_end)) {
                return xd;
            }
        } else if (index >= range_end) {
            if ((xd_end >= range_start) && (xd_end <= range_end) && (index == xd_start)) {
                return  xd;
            }
        } else  {
            if ((xd_end >= range_start) && (xd_end <= index) && (xd_start >= index) && (xd_start <= range_end)) {
                if (xd_end == index) {
                    return xd;
                }
                if (xd_start == index) {
                    return xd;
                }
                ++xd;
                continue;
            }
        }
        xd = compute_bigmin(xd, z_min, z_max);
    }
    return xd;
}