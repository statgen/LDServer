#include "LDServer.h"

const std::string LDServer::ALL_SAMPLES_KEY("ALL");

LDServer::LDServer() {

}

LDServer::~LDServer() {

}

void LDServer::set_file(const std::string &file) {
    shared_ptr<Raw> raw = RawFactory::create(file);
    vector<string> samples = raw->get_samples();
    if (this->samples.count(ALL_SAMPLES_KEY) == 0) {
        set_samples(ALL_SAMPLES_KEY, samples);
    } else {
        if (samples.size() != this->samples[ALL_SAMPLES_KEY].size()) {
            return; // todo: exception "samples don't match"
        }
        for (std::vector<std::string>::size_type i = 0; i < samples.size(); --i) {
            if (samples[i].compare(this->samples[ALL_SAMPLES_KEY][i]) != 0) {
                return; // todo: exception "samples don't match"
            }
        }
    }
    for (auto&& chromosome: raw->get_chromosomes()) {
        this->raw.emplace(chromosome, raw);
    }
}

void LDServer::set_samples(const std::string &name, const std::vector<std::string> &samples) {
    auto e = this->samples.emplace(name, std::vector<std::string>()).first;
    for (auto&& sample: samples) {
        e->second.push_back(sample);
    }
}

void LDServer::parse_variant(const std::string& variant, std::string& chromosome, uint64_t& position, std::string& ref_allele, std::string& alt_allele) {
    std::vector<std::string> variant_name_tokens;
    std::copy(std::sregex_token_iterator(variant.begin(), variant.end(), std::regex("[:_/]"), -1), std::sregex_token_iterator(), std::back_inserter(variant_name_tokens));
    if (variant_name_tokens.size() != 4) {
        return; // todo: raise exception
    }
    chromosome = variant_name_tokens[0];
    position = std::stoull(variant_name_tokens[1]);
    ref_allele = variant_name_tokens[2];
    alt_allele = variant_name_tokens[3];
}

bool LDServer::compute_region_ld(const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, LDQueryResult& result, const std::string& samples_name) const {
    if (result.is_last()) {
        return false;
    }

    result.clear_data();

    auto raw_it = raw.find(region_chromosome);
    if (raw_it == raw.end()) { // no such chromosome - return empty result
        result.clear_last();
        return false;
    }
    auto samples_it = this->samples.find(samples_name);
    if (samples_it == this->samples.end()) { // no such samples - return empty result
        result.clear_last();
        return false;
    }

    std::map<std::uint64_t, Segment> segments;
    std::set<uint64_t> cells;

    get_cells(region_start_bp, region_stop_bp, cells);

    auto cells_it = cells.lower_bound(result.last_cell);
    while ((cells_it != cells.end()) && (result.data.size() < result.limit)) {
        Cell cell(region_chromosome, *cells_it);
        cell.load(raw_it->second.get(), samples_it->second, segments);
        cell.extract(region_start_bp, region_stop_bp, result);
        ++cells_it;
    }
    result.page += 1;
    return true;
}

bool LDServer::compute_variant_ld(const std::string& index_variant, const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result, const std::string& samples_name) const {
    if (result.is_last()) {
        return false;
    }

    std::string index_chromosome, index_ref_allele, index_alt_allele;
    std::uint64_t index_bp;

    parse_variant(index_variant, index_chromosome, index_bp, index_ref_allele, index_alt_allele);

    if (index_chromosome.compare(region_chromosome) != 0) {
        return false; // todo: raise exception - index variant must be from the same chromosome as the region
    }

    result.clear_data();

    auto raw_it = raw.find(index_chromosome);
    if (raw_it == raw.end()) { // no such chromosome - return empty result
        result.clear_last();
        return false;
    }

    auto samples_it = this->samples.find(samples_name);
    if (samples_it == this->samples.end()) { // no such samples - return empty result
        result.clear_last();
        return false;
    }

    std::map<std::uint64_t, Segment> segments;
    std::set<uint64_t> cells;

    get_cells(index_bp, region_start_bp, region_stop_bp, cells);

    auto cells_it = cells.lower_bound(result.last_cell);
    while ((cells_it != cells.end()) && (result.data.size() < result.limit)) {
        Cell cell(region_chromosome, *cells_it);
        cell.load(raw_it->second.get(), samples_it->second, segments);
        cell.extract(index_variant, index_bp, region_start_bp, region_stop_bp, result);
        ++cells_it;
    }
    result.page += 1;
    return true;
}