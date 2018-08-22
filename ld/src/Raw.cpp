#include "Raw.h"

Raw::Raw(const string& file) : file(file) {

}

Raw::~Raw() {

}

RawVCF::~RawVCF() {
}

vector<string> RawVCF::get_samples() const {
    return savvy::vcf::reader<1>(file, savvy::fmt::gt).samples();
}

vector<string> RawVCF::get_chromosomes() const {
    return savvy::vcf::indexed_reader<1>(file, {""}, savvy::fmt::gt).chromosomes();
}

void RawVCF::load(const vector<string>& samples, Segment& segment) const {
    savvy::vcf::indexed_reader<1> f(file, {segment.chromosome, segment.start_bp, segment.stop_bp}, savvy::fmt::gt);
    f.subset_samples({samples.begin(), samples.end()});
    segment.sp_mat_rowind.clear();
    segment.sp_mat_colind.clear();
    segment.names.clear();
    segment.positions.clear();
    std::stringstream ss;
    savvy::site_info anno;
    savvy::armadillo::sparse_vector<float> alleles;
    while (f.read(anno, alleles)) {
        if (alleles.n_nonzero > 0) {
            ss.str("");
            ss << anno.chromosome() << ":" << anno.position() << "_" << anno.ref() << "/" << anno.alt();
            segment.names.emplace_back(ss.str());
            segment.positions.push_back(anno.position());
            segment.sp_mat_colind.push_back(segment.sp_mat_rowind.size());
            for (auto it = alleles.begin(); it != alleles.end(); ++it) {
                segment.sp_mat_rowind.push_back(it.row());
            }
        }
    }
    segment.sp_mat_colind.push_back(segment.sp_mat_rowind.size());
}

RawSAV::~RawSAV() {
}

vector<string> RawSAV::get_samples() const {
    return savvy::reader(file, savvy::fmt::gt).samples();
}

vector<string> RawSAV::get_chromosomes() const {
    return savvy::indexed_reader(file, {""}, savvy::fmt::gt).chromosomes();
}

void RawSAV::load(const vector<string>& samples, Segment& segment) const {
    savvy::indexed_reader f(file, {segment.chromosome, segment.start_bp, segment.stop_bp}, savvy::fmt::gt);
    f.subset_samples({samples.begin(), samples.end()});
    segment.sp_mat_rowind.clear();
    segment.sp_mat_colind.clear();
    segment.names.clear();
    segment.positions.clear();
    std::stringstream ss;
    savvy::site_info anno;
    savvy::armadillo::sparse_vector<float> alleles;
    while (f.read(anno, alleles)) {
        if (alleles.n_nonzero > 0) {
            ss.str("");
            segment.names.emplace_back(ss.str());
            segment.positions.push_back(anno.position());
            segment.sp_mat_colind.push_back(segment.sp_mat_rowind.size());
            for (auto it = alleles.begin(); it != alleles.end(); ++it) {
                segment.sp_mat_rowind.push_back(it.row());
            }
        }
    }
    segment.sp_mat_colind.push_back(segment.sp_mat_rowind.size());
}

shared_ptr<Raw> RawFactory::create(const string &file) {
    if ((file.length() >= 4) && (file.compare(file.length() - 4, 4, ".sav") == 0)) {
        return shared_ptr<RawSAV>(new RawSAV(file));
    } else if ((file.length() >= 7) && (file.compare(file.length() - 7, 7, ".vcf.gz") == 0)) {
        return shared_ptr<RawVCF>(new RawVCF(file));
    } else if ((file.length() >= 4) && (file.compare(file.length() - 4, 4, ".bcf") == 0)) {
        return shared_ptr<RawVCF>(new RawVCF(file));
    } else {
        //todo: throw exception "unrecognized file format"
    }
    return nullptr;
}

