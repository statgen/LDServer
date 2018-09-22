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
#include <savvy/reader.hpp>
#include <savvy/armadillo_vector.hpp>

using namespace std;

class Segment {
private:
    string key;
    bool cached;
    bool names_loaded;
    bool genotypes_loaded;


public:
    string chromosome;
    uint64_t start_bp;
    uint64_t stop_bp;
    uint64_t n_haplotypes;
    vector<string> names;
    vector<uint64_t> positions;

    vector<arma::uword> sp_mat_rowind;
    vector<arma::uword> sp_mat_colind;

    Segment(uint32_t unique_key, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp);
    Segment(Segment&& segment);
    virtual ~Segment();

    const char* get_key() const;
    uint64_t get_key_size() const;

    void clear();
    void clear_names();
    void clear_genotypes();
    void add(savvy::site_info& anno, savvy::armadillo::sparse_vector<float>& alleles);
    void add_name(savvy::site_info& anno, savvy::armadillo::sparse_vector<float>& alleles);
    void add_genotypes(savvy::armadillo::sparse_vector<float>& alleles);
    void freeze();
    void freeze_names();
    void freeze_genotypes();

    void load(redisContext* redis_cache);
    void save(redisContext* redis_cache);

    bool is_cached() const;
    bool has_names() const;
    bool has_genotypes() const;

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
