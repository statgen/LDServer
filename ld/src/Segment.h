#ifndef LDSERVER_REGION_H
#define LDSERVER_REGION_H

#include <iostream>
#include <vector>
#include <string>
#include <savvy/armadillo_vector.hpp>

using namespace std;

class Segment {
private:

public:
    string chromosome;
    uint64_t start_bp;
    uint64_t stop_bp;
    vector<arma::uword> sp_mat_rowind;
    vector<arma::uword> sp_mat_colind;
    vector<std::string> names;
    vector<std::uint64_t> positions;

    Segment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp);
    virtual ~Segment();
};


#endif
