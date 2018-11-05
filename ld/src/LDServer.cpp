#include "LDServer.h"

const std::string LDServer::ALL_SAMPLES_KEY("ALL");

LDServer::LDServer(uint32_t segment_size) : segment_size(segment_size), cache_enabled(false), cache_key(0u), cache_hostname(""), cache_port(0), cache_context(nullptr) {

}

LDServer::~LDServer() {
    if (cache_context != nullptr) {
        redisFree(cache_context);
        cache_context = nullptr;
    }
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

void LDServer::enable_cache(uint32_t cache_key, const string& hostname, int port) {
    if (cache_context == nullptr) {
        this->cache_key = cache_key;
        cache_hostname = hostname;
        cache_port = port;
        struct timeval timeout = {1, 500000}; // 1.5 seconds
        cache_context = redisConnectWithTimeout(cache_hostname.c_str(), cache_port, timeout);
        if (cache_context == nullptr) {
            return; // todo: exception "can't connect to Redis cache server" (see code in context->err)
        }
        cache_enabled = true;
    }
}

void LDServer::disable_cache() {
    if (cache_enabled) {
        if (cache_context != nullptr) {
            redisFree(cache_context);
            cache_context = nullptr;
        }
        cache_enabled = false;
    }
}

string LDServer::make_cell_cache_key(uint32_t cache_key, const string& samples_name, correlation correlation_type, const string& chromosome, uint64_t morton_code) {
    stringstream os(ios::binary | ios::out);
    os.write(reinterpret_cast<const char*>(&cache_key), sizeof(cache_key));
    os.write(samples_name.c_str(), samples_name.size());
    os.write(chromosome.c_str(), chromosome.size());
    os.write(reinterpret_cast<const char*>(&correlation_type), sizeof(correlation_type));
    os.write(reinterpret_cast<const char*>(&morton_code), sizeof(morton_code));
    os.flush();
    return os.str();
}

string LDServer::make_segment_cache_key(uint32_t cache_key, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp) {
    stringstream os(ios::binary | ios::out);
    os.write(reinterpret_cast<const char*>(&cache_key), sizeof(cache_key));
    os.write(samples_name.c_str(), samples_name.size());
    os.write(chromosome.c_str(), chromosome.size());
    os.write(reinterpret_cast<const char*>(&start_bp), sizeof(start_bp));
    os.write(reinterpret_cast<const char*>(&stop_bp), sizeof(stop_bp));
    os.flush();
    return os.str();
}

void LDServer::parse_variant(const std::string& variant, std::string& chromosome, uint64_t& position, std::string& ref_allele, std::string& alt_allele) {
    std::vector<std::string> variant_name_tokens;
    auto separator = std::regex("[:_/]");
    std::copy(std::sregex_token_iterator(variant.begin(), variant.end(), separator, -1), std::sregex_token_iterator(), std::back_inserter(variant_name_tokens));
    if (variant_name_tokens.size() != 4) {
        return; // todo: raise exception
    }
    chromosome = variant_name_tokens[0];
    position = std::stoull(variant_name_tokens[1]);
    ref_allele = variant_name_tokens[2];
    alt_allele = variant_name_tokens[3];
}

shared_ptr<Segment> LDServer::load_segment(const shared_ptr<Raw>& raw, const string& samples_name, bool only_variants, const std::string& chromosome, uint64_t i, std::map<std::uint64_t, shared_ptr<Segment>>& segments) const {
//    auto start = std::chrono::system_clock::now();
    uint64_t start_bp = i * segment_size;
    uint64_t stop_bp = i * segment_size + segment_size - 1u;
    string key = make_segment_cache_key(cache_key, samples_name, chromosome, start_bp, stop_bp);
    auto segment_it = segments.find(i);
    if (segment_it == segments.end()) {
        segment_it = segments.emplace(make_pair(i, make_shared<Segment>(chromosome, start_bp, stop_bp))).first;
        if (cache_enabled) {
            segment_it->second->load(cache_context, key);
        }
    }
    if (segment_it->second->is_cached()) {
        if ((!only_variants) && (!segment_it->second->has_genotypes())) {
            raw->load_genotypes(*(segment_it->second));
        }
    } else {
        if (only_variants) {
            if (!segment_it->second->has_names()) {
                raw->load_names(*(segment_it->second));
            }
        } else {
            if (!segment_it->second->has_names() && !segment_it->second->has_genotypes()) {
                raw->load(*(segment_it->second));
            } else if (!segment_it->second->has_names()) {
                raw->load_names(*(segment_it->second));
            } else if (!segment_it->second->has_genotypes()){
                raw->load_genotypes(*(segment_it->second));
            }
        }
        if (cache_enabled) {
            segment_it->second->save(cache_context, key);
        }
    }
//    auto end = std::chrono::system_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    std::cout << "Segment load elapsed time: " << elapsed.count() << " s\n";
    return segment_it->second;
}

vector<string> LDServer::get_chromosomes() {
    vector<string> chromosomes;
    for (const auto& x : this->raw) {
        chromosomes.emplace_back(x.first);
    }
    return chromosomes;
}

bool LDServer::compute_region_ld(const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, correlation correlation_type, LDQueryResult& result, const std::string& samples_name) const {
//    auto start = std::chrono::system_clock::now();
    if (result.is_last()) {
        return false;
    }

    result.clear_data();

    result.page += 1;

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

    raw_it->second->open(region_chromosome, samples_it->second);

    std::map<std::uint64_t, shared_ptr<Segment>> segments;

    uint64_t segment_i = region_start_bp / segment_size;
    uint64_t segment_j = region_stop_bp / segment_size;
    uint64_t z_min = to_morton_code(segment_i, segment_i);
    uint64_t z_max = to_morton_code(segment_j, segment_j);

    uint64_t z = get_next_z(segment_i, segment_j, z_min, z_max, result.last_cell > z_min ? result.last_cell : z_min);
    uint64_t i = 0u, j = 0u;

    CellFactory factory;

    while (z <= z_max) {
        string key = make_cell_cache_key(cache_key, samples_name, correlation_type, region_chromosome, z);
        from_morton_code(z, i, j);
        shared_ptr<Cell> cell = factory.create(correlation_type, i, j);
        if (cache_enabled) {
            cell->load(cache_context, key);
        }
        cell->segment_i = load_segment(raw_it->second, samples_name, cell->is_cached(), region_chromosome, cell->get_i(), segments);
        if (!cell->is_diagonal()) {
            cell->segment_j = load_segment(raw_it->second, samples_name, cell->is_cached(), region_chromosome, cell->get_j(), segments);
        }
        if (!cell->is_cached()) {
            cell->compute();
            if (cache_enabled) {
                cell->save(cache_context, key);
            }
        }
        cell->extract(region_start_bp, region_stop_bp, result);
        if ((result.last_i >= 0) && (result.last_j >= 0)) {
            result.last_cell = z;
            break;
        }
        z = get_next_z(segment_i, segment_j, z_min, z_max, ++z);
        if (result.data.size() >= result.limit) {
            if (z <= z_max) {
                result.last_cell = z;
                result.last_i = 0;
                result.last_j = 0;
            }
            break;
        }
    }
//    auto end = std::chrono::system_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    std::cout << "Page elapsed time: " << elapsed.count() << " s\n";
    return true;
}

bool LDServer::compute_variant_ld(const std::string& index_variant, const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, correlation correlation_type, struct LDQueryResult& result, const std::string& samples_name) const {
    if (result.is_last()) {
        return false;
    }

    result.page += 1;

    std::string index_chromosome, index_ref_allele, index_alt_allele;
    std::uint64_t index_bp;

    parse_variant(index_variant, index_chromosome, index_bp, index_ref_allele, index_alt_allele);

    if (index_chromosome.compare(region_chromosome) != 0) {
        result.clear_last();
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

    raw_it->second->open(index_chromosome, samples_it->second);

    std::map<std::uint64_t, shared_ptr<Segment>> segments;

    uint64_t segment_index = index_bp / segment_size;
    uint64_t segment_i = region_start_bp / segment_size;
    uint64_t segment_j = region_stop_bp / segment_size;
    uint64_t z_min = 0,z_max = 0;
    if (segment_index <= segment_i) { // only upper triangle of the matrix must be used
        z_min = to_morton_code(segment_i, segment_index);
        z_max = to_morton_code(segment_j, segment_index);
    } else if ((segment_index > segment_i) && (segment_index < segment_j)) {
        z_min = to_morton_code(segment_index, segment_i);
        z_max = to_morton_code(segment_j, segment_index);
    } else {
        z_min = to_morton_code(segment_index, segment_i);
        z_max = to_morton_code(segment_index, segment_j);
    }
    uint64_t z = get_next_z(segment_index, segment_i, segment_j, z_min, z_max, result.last_cell > z_min ? result.last_cell : z_min);
    uint64_t i = 0u, j = 0u;
    CellFactory factory;
    while (z <= z_max) {
        string key = make_cell_cache_key(cache_key, samples_name, correlation_type, region_chromosome, z);
        from_morton_code(z, i, j);
        shared_ptr<Cell> cell = factory.create(correlation_type, i, j);
        if (cache_enabled) {
            cell->load(cache_context, key);
        }
        cell->segment_i = load_segment(raw_it->second, samples_name, cell->is_cached(), region_chromosome, cell->get_i(), segments);
        if (!cell->is_diagonal()) {
            cell->segment_j = load_segment(raw_it->second, samples_name, cell->is_cached(), region_chromosome, cell->get_j(), segments);
        }
        if (!cell->is_cached()) {
            cell->compute();
            if (cache_enabled) {
                cell->save(cache_context, key);
            }
        }
        cell->extract(index_variant, index_bp, region_start_bp, region_stop_bp, result);
        if (result.last_j >= 0) {
            result.last_cell = z;
            break;
        }
        z = get_next_z(segment_index, segment_i, segment_j, z_min, z_max, ++z);
        if (result.data.size() >= result.limit) {
            if (z <= z_max) {
                result.last_cell = z;
                result.last_j = 0;
            }
            break;
        }
    }
    return true;
}