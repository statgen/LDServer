#ifndef LDSERVER_REGION_H
#define LDSERVER_REGION_H

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
    char* key;
    uint64_t key_size;

    bool has_genotypes;

public:
    string chromosome;
    uint64_t start_bp;
    uint64_t stop_bp;
    vector<string> names;
    vector<uint64_t> positions;

    vector<arma::uword> sp_mat_rowind;
    vector<arma::uword> sp_mat_colind;

    Segment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp);
    Segment(Segment&& segment);
    virtual ~Segment();

    const char* get_key() const;
    uint64_t get_key_size() const;

    template <class Archive>
    void save( Archive & ar ) const
    {
        ar( names, positions );
    }

    template <class Archive>
    void load( Archive & ar )
    {
        ar( names, positions );
    }
};

#endif
