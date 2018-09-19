#ifndef LDSERVER_SEGMENT_H
#define LDSERVER_SEGMENT_H

#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <strstream>
#include <savvy/armadillo_vector.hpp>
#include <hiredis/hiredis.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

using namespace std;

class Segment {
private:
    string key;
    bool cached;

public:
    string chromosome;
    uint64_t start_bp;
    uint64_t stop_bp;
    uint64_t n_haplotypes;
    vector<string> names;
    vector<uint64_t> positions;

    vector<arma::uword> sp_mat_rowind;
    vector<arma::uword> sp_mat_colind;

    bool genotypes_loaded;
    bool variants_loaded;

    Segment(uint32_t unique_key, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp);
    Segment(Segment&& segment);
    virtual ~Segment();

    const char* get_key() const;
    uint64_t get_key_size() const;

    void load(redisContext* redis_cache);
    void save(redisContext* redis_cache);

    bool is_cached() const;

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

#endif
