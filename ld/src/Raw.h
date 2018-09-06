#ifndef LDSERVER_RAW_H
#define LDSERVER_RAW_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <savvy/reader.hpp>
#include <savvy/armadillo_vector.hpp>
#include "Segment.h"

using namespace std;

class Raw {
protected:
    string file;

public:
    Raw(const string& file);
    virtual ~Raw();

    virtual vector<string> get_samples() const = 0;
    virtual vector<string> get_chromosomes() const = 0;
    virtual void load(const vector<string>& samples, const shared_ptr<Segment>& segment) const = 0;
};

class RawVCF : public Raw {
public:
    using Raw::Raw;
    virtual ~RawVCF();

    virtual vector<string> get_samples() const;
    virtual vector<string> get_chromosomes() const;
    virtual void load(const vector<string>& samples, const shared_ptr<Segment>& segment) const;
};

class RawSAV : public Raw {
public:
    using Raw::Raw;
    virtual ~RawSAV();

    virtual vector<string> get_samples() const;
    virtual vector<string> get_chromosomes() const;
    virtual void load(const vector<string>& samples, const shared_ptr<Segment>& segment) const;
};

class RawFactory {
public:
    static shared_ptr<Raw> create(const string& file);
};


#endif
