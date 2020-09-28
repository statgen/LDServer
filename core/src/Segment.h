#ifndef LDSERVER_SEGMENT_H
#define LDSERVER_SEGMENT_H

#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <savvy/armadillo_vector.hpp>
#include <hiredis/hiredis.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <savvy/reader.hpp>
#include <armadillo>
#include "Types.h"

using namespace std;

enum genotypes_store : uint8_t {
    CSC_ALL_ONES, // all-ones matrix in compressed sparse column representation,
    CSC, // matrix in compressed sparse column representation
    BITSET // bitsets
};

class Segment {
protected:
    bool cached;
    bool names_loaded;
    bool genotypes_loaded;

    string chromosome;
    uint64_t start_bp;
    uint64_t stop_bp;
    uint64_t n_haplotypes;

    // Vector of variant IDs in EPACTS format (chr:pos_ref/alt).
    vector<string> names;
    vector<uint64_t> positions;

    // tells in which format to store the genotypes
    genotypes_store store;

    // for CSC representation
    vector<arma::uword> sp_mat_rowind;
    vector<arma::uword> sp_mat_colind;
    vector<float> sp_mat_values;
    vector<double> means;
    bool nans = false;

    // for BITSET representation
    vector<float> freqs;
    vector<vector<bool>> alleles;
    vector<vector<unsigned int>> alt_carriers;

public:
    Segment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp, genotypes_store store);
    Segment(Segment&& segment);
    virtual ~Segment();

    void clear();
    void clear_names();
    void clear_genotypes();
    void add(savvy::site_info& anno, savvy::compressed_vector<float>& alleles);
    void add_name(savvy::site_info& anno, savvy::compressed_vector<float>& alleles);
    void add_genotypes(savvy::compressed_vector<float>& alleles);
    void freeze();
    void freeze_names();
    void freeze_genotypes();

    bool is_empty() const;
    bool is_cached() const;
    bool has_names() const;
    bool has_genotypes() const;

    const char* get_key() const;
    uint64_t get_key_size() const;
    const string& get_chromosome() const;
    uint64_t get_start_bp() const;
    uint64_t get_stop_bp() const;
    uint64_t get_n_haplotypes() const;
    uint64_t get_n_genotypes() const;
    uint32_t get_n_variants() const;
    const string& get_name(int i) const;
    uint64_t get_position(int i) const;
    arma::sp_fmat get_genotypes();
    const vector<float>& get_freqs() const;
    const vector<vector<bool>>& get_alleles() const;
    const vector<vector<unsigned int>>& get_alt_carriers() const;
    genotypes_store get_store() const;

    static void create_pair(Segment& segment1, Segment& segment2, int index_i, int index_j, double value, LDQueryResult& result);
    bool overlaps_region(uint64_t region_start_bp, uint64_t region_stop_bp, int& from_index, int& to_index) const;
    bool overlaps_variant(const string& name, uint64_t bp, int& index) const;

    inline bool has_nans() const { return nans; }
    inline const vector<double>& get_means() const { return means; }

    /**
     * Load/save functions for redis.
     * These are also overloaded below for loading/saving from serialized binary.
     * @param redis_cache
     * @param key
     */
    virtual void load(redisContext* redis_cache, const string& key);
    virtual void save(redisContext* redis_cache, const string& key);

    /**
     * Load/save functions for binary format.
     * Only the number of haplotypes, the list of variant EPACTS IDs, and their positions on the chromosome are stored.
     * @tparam Archive
     * @param ar
     */
    template <class Archive>
    void load( Archive & ar )
    {
        ar( n_haplotypes, names, positions );
    }

    template <class Archive>
    void save( Archive & ar ) const
    {
        ar( n_haplotypes, names, positions );
    }
};

typedef shared_ptr<vector<shared_ptr<Segment>>> SharedSegmentVector;
inline SharedSegmentVector make_shared_segment_vector() { return make_shared<vector<shared_ptr<Segment>>>(); }

#endif
