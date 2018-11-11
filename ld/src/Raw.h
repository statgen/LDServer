#ifndef LDSERVER_RAW_H
#define LDSERVER_RAW_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <savvy/reader.hpp>
#include "Segment.h"

using namespace std;

class Raw {
protected:
    string file;

public:
    Raw(const string& file);
    virtual ~Raw();

    virtual void open(const string& chromosome, const vector<string>& samples, bool coded012) = 0;
    virtual vector<string> get_samples() const = 0;
    virtual vector<string> get_chromosomes() const = 0;
    virtual void load(Segment& segment) = 0;
    virtual void load_names(Segment &segment) = 0;
    virtual void load_genotypes(Segment &segment) = 0;
};

class RawVCF : public Raw {
private:
    unique_ptr<savvy::vcf::indexed_reader<1>> f;
    bool has_cached;
    savvy::site_info anno;
    savvy::compressed_vector<float> alleles;

public:
    using Raw::Raw;
    virtual ~RawVCF();

    void open(const string& chromosome, const vector<string>& samples, bool coded012) override;
    vector<string> get_samples() const override;
    vector<string> get_chromosomes() const override;
    void load(Segment& segment) override;
    void load_names(Segment &segment) override;
    void load_genotypes(Segment &segment) override;
};

class RawSAV : public Raw {
private:
    unique_ptr<savvy::indexed_reader> f;
    bool has_cached;
    savvy::site_info anno;
    savvy::compressed_vector<float> alleles;

public:
    using Raw::Raw;
    virtual ~RawSAV();

    void open(const string& chromosome, const vector<string>& samples, bool coded012) override;
    vector<string> get_samples() const override;
    vector<string> get_chromosomes() const override;
    void load(Segment& segment) override;
    void load_names(Segment &segment) override;
    void load_genotypes(Segment &segment) override;
};

class RawFactory {
public:
    static shared_ptr<Raw> create(const string& file);
};


#endif
