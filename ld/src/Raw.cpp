#include "Raw.h"

Raw::Raw(const string& file) : file(file) {

}

Raw::~Raw() {

}

RawVCF::~RawVCF() {

}

void RawVCF::open(const string& chromosome, const vector<string>& samples, bool coded012) {
    f.release();
    f = unique_ptr<savvy::vcf::indexed_reader<1>>(new savvy::vcf::indexed_reader<1>(file, {chromosome}, coded012 ? savvy::fmt::ac : savvy::fmt::gt));
    f->subset_samples({samples.begin(), samples.end()});
    has_cached = false;
}

vector<string> RawVCF::get_samples() const {
    return savvy::vcf::reader<1>(file, savvy::fmt::gt).samples();
}

vector<string> RawVCF::get_chromosomes() const {
    return savvy::vcf::indexed_reader<1>(file, {""}, savvy::fmt::gt).chromosomes();
}

void RawVCF::load(Segment& segment) {
    segment.clear();
    if (has_cached && (segment.get_start_bp() <= anno.position()) && (anno.position() <= segment.get_stop_bp())) {
        segment.add(anno, alleles);
    } else{
        f->reset_region({segment.get_chromosome(), segment.get_start_bp(), numeric_limits<int>::max() - 1});
    }
    has_cached = false;
    while (f->read(anno, alleles).good()) {
        if (anno.position() > segment.get_stop_bp()) {
            has_cached = true;
            break;
        }
        segment.add(anno, alleles);
    }
    segment.freeze();
}

void RawVCF::load_names(Segment &segment) {
    segment.clear_names();
    if (has_cached && (segment.get_start_bp() <= anno.position()) && (anno.position() <= segment.get_stop_bp())) {
        segment.add_name(anno, alleles);
    } else{
        f->reset_region({segment.get_chromosome(), segment.get_start_bp(), numeric_limits<int>::max() - 1});
    }
    has_cached = false;
    while (f->read(anno, alleles)) {
        if (anno.position() > segment.get_stop_bp()) {
            has_cached = true;
            break;
        }
        segment.add_name(anno, alleles);
    }
    segment.freeze_names();
}

void RawVCF::load_genotypes(Segment &segment) {
    segment.clear_genotypes();
    if (has_cached && (segment.get_start_bp() <= anno.position()) && (anno.position() <= segment.get_stop_bp())) {
        segment.add_genotypes(alleles);
    } else{
        f->reset_region({segment.get_chromosome(), segment.get_start_bp(), numeric_limits<int>::max() - 1});
    }
    while (f->read(anno, alleles)) {
        if (anno.position() > segment.get_stop_bp()) {
            has_cached = true;
            break;
        }
        segment.add_genotypes(alleles);
    }
    segment.freeze_genotypes();
}

RawSAV::~RawSAV() {
}

void RawSAV::open(const string& chromosome, const vector<string>& samples, bool coded012) {
    f.release();
    f = unique_ptr<savvy::indexed_reader>(new savvy::indexed_reader(file, {chromosome}, coded012 ? savvy::fmt::ac : savvy::fmt::gt));
    f->subset_samples({samples.begin(), samples.end()});
    has_cached = false;
}

vector<string> RawSAV::get_samples() const {
    return savvy::reader(file, savvy::fmt::gt).samples();
}

vector<string> RawSAV::get_chromosomes() const {
    return savvy::indexed_reader(file, {""}, savvy::fmt::gt).chromosomes();
}

void RawSAV::load(Segment& segment) {
    segment.clear();
    if (has_cached && (segment.get_start_bp() <= anno.position()) && (anno.position() <= segment.get_stop_bp())) {
        segment.add(anno, alleles);
    } else {
        f->reset_region({segment.get_chromosome(), segment.get_start_bp()});
    }
    has_cached = false;
    while (f->read(anno, alleles).good()) {
        if (anno.position() > segment.get_stop_bp()) {
            has_cached = true;
            break;
        }
        segment.add(anno, alleles);
    }
    segment.freeze();
}

void RawSAV::load_names(Segment &segment) {
    segment.clear_names();
    if (has_cached && (segment.get_start_bp() <= anno.position()) && (anno.position() <= segment.get_stop_bp())) {
        segment.add_name(anno, alleles);
    } else{
        f->reset_region({segment.get_chromosome(), segment.get_start_bp()});
    }
    has_cached = false;
    while (f->read(anno, alleles)) {
        if (anno.position() > segment.get_stop_bp()) {
            has_cached = true;
            break;
        }
        segment.add_name(anno, alleles);
    }
    segment.freeze_names();
}

void RawSAV::load_genotypes(Segment &segment) {
    segment.clear_genotypes();
    if (has_cached && (segment.get_start_bp() <= anno.position()) && (anno.position() <= segment.get_stop_bp())) {
        segment.add_genotypes(alleles);
    } else{
        f->reset_region({segment.get_chromosome(), segment.get_start_bp()});
    }
    while (f->read(anno, alleles)) {
        if (anno.position() > segment.get_stop_bp()) {
            has_cached = true;
            break;
        }
        segment.add_genotypes(alleles);
    }
    segment.freeze_genotypes();
}

shared_ptr<Raw> RawFactory::create(const string &file) {
    if ((file.length() >= 4) && (file.compare(file.length() - 4, 4, ".sav") == 0)) {
        return shared_ptr<RawSAV>(new RawSAV(file));
    } else if ((file.length() >= 7) && (file.compare(file.length() - 7, 7, ".vcf.gz") == 0)) {
        return shared_ptr<RawVCF>(new RawVCF(file));
    } else if ((file.length() >= 4) && (file.compare(file.length() - 4, 4, ".bcf") == 0)) {
        return shared_ptr<RawVCF>(new RawVCF(file));
    }
    throw runtime_error("Unknown genotype file type");
}

