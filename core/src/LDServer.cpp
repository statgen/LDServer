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

void LDServer::force_samples(const std::string &name, const std::vector<std::string> &samples) {
    this->samples[name] = samples;
}

/**
 * Enable the cache
 * @param cache_key - This key gets hashed together with a number of other parameters to create segment and cell cache keys.
 *                    It is usually the reference panel (or genotype dataset) ID.
 * @param hostname
 * @param port
 */
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

/**
 * Load a segment of genotypes.
 * @param raw
 * @param store
 * @param samples_name
 * @param only_variants This corresponds to cell->is_cached(). It corresponds to "do we need to load only the variant names?"
 * @param chromosome
 * @param i
 * @param segments
 * @return
 */
shared_ptr<Segment> LDServer::load_segment(const shared_ptr<Raw>& raw, genotypes_store store, const string& samples_name, bool only_variants, const std::string& chromosome, uint64_t i, std::map<std::uint64_t, shared_ptr<Segment>>& segments) const {
//    auto start = std::chrono::system_clock::now();
    uint64_t start_bp = i * segment_size;
    uint64_t stop_bp = i * segment_size + segment_size - 1u;
    string key = make_segment_cache_key(cache_key, samples_name, chromosome, start_bp, stop_bp);
    auto segment_it = segments.find(i);
    if (segment_it == segments.end()) {
        segment_it = segments.emplace(make_pair(i, make_shared<Segment>(chromosome, start_bp, stop_bp, store))).first;
        if (cache_enabled) {
            segment_it->second->load(cache_context, key);
        }
    }
    if (segment_it->second->is_cached()) {
        // If the segment was cached, but the cell was not cached, we need to load the genotypes in order to perform
        // the correlation calculations. When a segment is cached, only the variant names/positions are stored, not
        // the genotypes.
        if ((!only_variants) && (!segment_it->second->has_genotypes())) {
            raw->load_genotypes(*(segment_it->second));
        }
    } else {
        // We only need to load the variant names, because this cell (and the computed corr values) were
        // previously cached, and so it isn't necessary to load the genotypes (no calculations will be performed.)
        if (only_variants) {
            if (!segment_it->second->has_names()) {
                raw->load_names(*(segment_it->second));
            }
        } else {
            // In this branch, the cell wasn't cached, and the segment was not cached either.
            // We may need to load both variant IDs/names and genotypes.
            if (!segment_it->second->has_names() && !segment_it->second->has_genotypes()) {
                // This loads both variant IDs / positions and also the genotypes.
                raw->load(*(segment_it->second));
            } else if (!segment_it->second->has_names()) {
                // Only names are loaded here, we already have the genotypes loaded.
                raw->load_names(*(segment_it->second));
            } else if (!segment_it->second->has_genotypes()){
                // Names were loaded from cache, we only need to load genotypes.
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

void LDServer::add_overlap_position(const uint64_t& pos) {
  uint64_t seg_i = pos / segment_size;
  allowed_segments.insert(seg_i);
}

template<template <typename...> class C, typename T>
bool contains(const C<T>& container, const T& v) {
  return container.find(v) != container.end();
}

/**
 * Return a vector of 2-tuples for the cartesian product of a container with itself.
 * For example, a container with elements A, B, C would result in: (A, A), (A, B), (A, C), (B, B), (B, C), (C, C)
 * @tparam C container type
 * @tparam T element type
 * @param container of elements of type T
 * @return shared pointer to a vector of tuples<T, T>, one tuple per combination
 */
template<template <typename...> class C, typename T>
shared_ptr<vector<tuple<T, T>>> product(const C<T>& container) {
  auto result = make_shared<vector<tuple<T, T>>>();
  for (auto it1 = container.begin(); it1 != container.end(); it1++) {
    auto it2(it1);
    for (; it2 != container.end(); it2++) {
      auto t = make_tuple(*it1, *it2);
      result->emplace_back(t);
    }
  }
  return result;
}

/**
 * Function to compute LD within a region.
 * @param region_chromosome
 * @param region_start_bp
 * @param region_stop_bp
 * @param correlation_type
 * @param result
 * @param samples_name
 * @param segments_out - Pass out a list of segments if a shared_ptr to a vector of segments was provided. This allows
 *   passing on segments to another function to be used for computations.
 * @return
 */
bool LDServer::compute_region_ld(const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, correlation correlation_type, LDQueryResult& result, const std::string& samples_name, bool diagonal, shared_ptr<vector<shared_ptr<Segment>>> segments_out) const {
//    auto start = std::chrono::system_clock::now();
    if (result.is_last()) {
        // It seems as though we reach this code if:
        //   1. Extraction of the final cell completed (last_i and last_j are -1)
        //   2. The while z <= z_max loop finished
        //   3. We previously returned at least 1 page of data
        return false;
    }

    result.clear_data();

    result.page += 1;

    auto raw_it = raw.find(region_chromosome);
    if (raw_it == raw.end()) { // no such chromosome - return empty result
        result.clear_last();
        std::cout << "Couldn't find requested chromosome: " + region_chromosome << std::endl;
        return false;
    }
    auto samples_it = this->samples.find(samples_name);
    if (samples_it == this->samples.end()) { // no such samples - return empty result
        result.clear_last();
        std::cout << "Couldn't find sample subset named: " + samples_name << std::endl;
        return false;
    }

    raw_it->second->open(region_chromosome, samples_it->second, correlation_type == correlation::COV);
    genotypes_store store = genotypes_store::CSC_ALL_ONES;
    if (correlation_type == correlation::COV) {
        store = genotypes_store::CSC;
    } else if (correlation_type == correlation::LD_RSQUARE_APPROX) {
        store = genotypes_store::BITSET;
    }

    std::map<std::uint64_t, shared_ptr<Segment>> segments;

    uint64_t segment_i = region_start_bp / segment_size;
    uint64_t segment_j = region_stop_bp / segment_size;
    uint64_t z_min = to_morton_code(segment_i, segment_i);
    uint64_t z_max = to_morton_code(segment_j, segment_j);

    uint64_t z = get_next_z(segment_i, segment_j, z_min, z_max, result.last_cell > z_min ? result.last_cell : z_min);
    uint64_t i = 0u, j = 0u;

    CellFactory factory;

    // If there's only a few segments to be used for computation, we only need to evaluate cells at morton codes
    // for each combination of those segments.
    set<uint64_t> allowed_z;
    if (!allowed_segments.empty()) {
      auto segment_combos = product(allowed_segments);
      for (auto&& pair : *segment_combos) {
        uint64_t zi = to_morton_code(get<0>(pair), get<1>(pair));
        if (zi >= z) {
          allowed_z.emplace(zi);
        }
      }
    }

    auto zit = allowed_z.begin();
    while (z <= z_max) {
        // This converts from the linear index (morton code) z into coordinates i and j (and stores to them directly.)
        from_morton_code(z, i, j);

        string key = make_cell_cache_key(cache_key, samples_name, correlation_type, region_chromosome, z);

        shared_ptr<Cell> cell = factory.create(correlation_type, i, j);
        if (cache_enabled) {
            cell->load(cache_context, key);
        }
        cell->segment_i = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_i(), segments);
        if (!cell->is_diagonal()) {
            cell->segment_j = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_j(), segments);
        }
        if (!cell->is_cached()) {
            cell->compute();
            if (cache_enabled) {
                cell->save(cache_context, key);
            }
        }

        cell->extract(region_start_bp, region_stop_bp, result, diagonal);

        // Setting last_i and last_j to -1 is the sentinel for flagging that extraction of this cell is either complete,
        // or had to stop due to the result object becoming full.
        if ((result.last_i >= 0) && (result.last_j >= 0)) {
            // We will need to reload this cell z next time, since we stopped in the middle of the segments somewhere.
            result.last_cell = z;
            break;
        }

        // We didn't fill up the result object yet, so we can continue to another cell.
        if (!allowed_z.empty()) {
          // We're iterating over a known series of morton codes, probably because we only wanted to calculate
          // over certain segments (like aggregation test covariance matrices.)
          zit++;
          if (zit != allowed_z.end()) {
            // There are more morton codes to iterate over
            z = *zit;
          }
          else {
            // We're done.
            z = z_max + 1;
          }
        }
        else {
          // This function finds the next morton code to use, with the knowledge that we only need the upper triangle
          // of the matrix.
          z = get_next_z(segment_i, segment_j, z_min, z_max, ++z);
        }

        if (result.n_correlations >= result.limit) {
            // If we reach this point, the following is true:
            //   1. The result object is full
            //   2. Cell extraction stopped because the cell had been parsed completely
            if (z <= z_max) {
                // We have more cells to parse, but the current one is completely finished. Set the resume point
                // to (0,0) of the next cell. Remember we "incremented" z above with get_next_z().
                result.last_cell = z;
                result.last_i = 0;
                result.last_j = 0;
            }
            break;
        }
    }

//  Before returning the result, fill out full variant names, chromosomes, and positions
    if (result.raw_variants.size() > 0) {
        result.data.variants.reserve(result.raw_variants.size());
        result.data.chromosomes.reserve(result.raw_variants.size());
        result.data.positions.reserve(result.raw_variants.size());
        for (auto&& entry: result.raw_variants) {
            auto segment = segments[entry.first];
            result.data.variants.push_back(segment->get_name(entry.second));
            result.data.chromosomes.push_back(segment->get_chromosome());
            result.data.positions.push_back(segment->get_position(entry.second));
        }
    }

    // If a pointer to a container of segments was provided, write out the segments we loaded here to it.
    // This can be used by the ScoreServer to calculate score statistics without reloading the genotypes.
    if (segments_out != nullptr) {
        for (auto&& kv : segments) {
            segments_out->emplace_back(kv.second);
        }
    }

//    auto end = std::chrono::system_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//    std::cout << "Page elapsed time: " << elapsed.count() << " s\n";
    return true;
}

//bool LDServer::compute_variant_ld(const std::string& index_variant, const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, correlation correlation_type, struct SingleVariantLDQueryResult& result, const std::string& samples_name) const {
//    if (result.is_last()) {
//        return false;
//    }
//
//    result.page += 1;
//
//    std::string index_chromosome, index_ref_allele, index_alt_allele;
//    std::uint64_t index_bp;
//
//    parse_variant(index_variant, index_chromosome, index_bp, index_ref_allele, index_alt_allele);
//
//    if (index_chromosome.compare(region_chromosome) != 0) {
//        result.clear_last();
//        return false; // todo: raise exception - index variant must be from the same chromosome as the region
//    }
//
//    result.clear_data();
//
//    auto raw_it = raw.find(index_chromosome);
//    if (raw_it == raw.end()) { // no such chromosome - return empty result
//        result.clear_last();
//        return false;
//    }
//
//    auto samples_it = this->samples.find(samples_name);
//    if (samples_it == this->samples.end()) { // no such samples - return empty result
//        result.clear_last();
//        return false;
//    }
//
//    raw_it->second->open(index_chromosome, samples_it->second, correlation_type == correlation::COV);
//    genotypes_store store = genotypes_store::CSC_ALL_ONES;
//    if (correlation_type == correlation::COV) {
//        store = genotypes_store::CSC;
//    } else if (correlation_type == correlation::LD_RSQUARE_APPROX) {
//        store = genotypes_store::BITSET;
//    }
//
//    std::map<std::uint64_t, shared_ptr<Segment>> segments;
//
//    uint64_t segment_index = index_bp / segment_size;
//    uint64_t segment_i = region_start_bp / segment_size;
//    uint64_t segment_j = region_stop_bp / segment_size;
//    uint64_t z_min = 0,z_max = 0;
//    if (segment_index <= segment_i) { // only upper triangle of the matrix must be used
//        z_min = to_morton_code(segment_i, segment_index);
//        z_max = to_morton_code(segment_j, segment_index);
//    } else if ((segment_index > segment_i) && (segment_index < segment_j)) {
//        z_min = to_morton_code(segment_index, segment_i);
//        z_max = to_morton_code(segment_j, segment_index);
//    } else {
//        z_min = to_morton_code(segment_index, segment_i);
//        z_max = to_morton_code(segment_index, segment_j);
//    }
//    uint64_t z = get_next_z(segment_index, segment_i, segment_j, z_min, z_max, result.last_cell > z_min ? result.last_cell : z_min);
//    uint64_t i = 0u, j = 0u;
//    CellFactory factory;
//    while (z <= z_max) {
//        string key = make_cell_cache_key(cache_key, samples_name, correlation_type, region_chromosome, z);
//        from_morton_code(z, i, j);
//        shared_ptr<Cell> cell = factory.create(correlation_type, i, j);
//        if (cache_enabled) {
//            cell->load(cache_context, key);
//        }
//        cell->segment_i = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_i(), segments);
//        if (!cell->is_diagonal()) {
//            cell->segment_j = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_j(), segments);
//        }
//        if (!cell->is_cached()) {
//            cell->compute();
//            if (cache_enabled) {
//                cell->save(cache_context, key);
//            }
//        }
//        cell->extract(index_variant, index_bp, region_start_bp, region_stop_bp, result);
//        if (result.last_j >= 0) {
//            result.last_cell = z;
//            break;
//        }
//        z = get_next_z(segment_index, segment_i, segment_j, z_min, z_max, ++z);
//        if (result.n_correlations >= result.limit) {
//            if (z <= z_max) {
//                result.last_cell = z;
//                result.last_j = 0;
//            }
//            break;
//        }
//    }
//
//    //  Before returning the result, fill out full variant names, chromosomes, and positions
//    if (result.raw_index_variant.second >= 0) {
//        auto segment = segments[result.raw_index_variant.first];
//        result.data.index_variant = segment->get_name(result.raw_index_variant.second);
//        result.data.index_chromosome = segment->get_chromosome();
//        result.data.index_position = segment->get_position(result.raw_index_variant.second);
//    }
//
//    if (result.raw_variants.size() > 0) {
//        result.data.variants.reserve(result.raw_variants.size());
//        result.data.chromosomes.reserve(result.raw_variants.size());
//        result.data.positions.reserve(result.raw_variants.size());
//        for (auto&& entry: result.raw_variants) {
//            auto segment = segments[entry.first];
//            result.data.variants.push_back(segment->get_name(entry.second));
//            result.data.chromosomes.push_back(segment->get_chromosome());
//            result.data.positions.push_back(segment->get_position(entry.second));
//        }
//    }
//
//    return true;
//}

bool LDServer::compute_variant_ld(const std::string& index_variant, const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, correlation correlation_type, struct SingleVariantLDQueryResult& result, const std::string& samples_name) const {
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

    raw_it->second->open(index_chromosome, samples_it->second, correlation_type == correlation::COV);
    genotypes_store store = genotypes_store::CSC_ALL_ONES;
    if (correlation_type == correlation::COV) {
        store = genotypes_store::CSC;
    } else if (correlation_type == correlation::LD_RSQUARE_APPROX) {
        store = genotypes_store::BITSET;
    }

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

//    uint64_t z = get_next_z(segment_index, segment_i, segment_j, z_min, z_max, result.last_cell > z_min ? result.last_cell : z_min);
    uint64_t i = 0u, j = 0u;
    CellFactory factory;

    uint64_t z_init = result.last_cell > z_min ? result.last_cell : z_min;
    unsigned int max_n_lookahead = 16u;
    if (omp_get_max_threads() < max_n_lookahead) {
        max_n_lookahead = omp_get_max_threads();
    }
    vector<tuple<uint64_t, string, shared_ptr<Cell>>> cell_fifo;
    bool filled = false;
    while (!filled) {
        cell_fifo.clear();

        for (unsigned int i_lookahead = 0u; i_lookahead < max_n_lookahead; ++i_lookahead) {
            uint64_t z = get_next_z(segment_index, segment_i, segment_j, z_min, z_max, z_init);
            string key = make_cell_cache_key(cache_key, samples_name, correlation_type, region_chromosome, z);
            from_morton_code(z, i, j);
            shared_ptr<Cell> cell = factory.create(correlation_type, i, j);
            if (cache_enabled) {
                cell->load(cache_context, key);
            }
            cell->segment_i = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_i(), segments);
            if (!cell->is_diagonal()) {
                cell->segment_j = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_j(), segments);
            }
            cell_fifo.emplace_back(z, move(key), move(cell));
            if (z == z_max) {
                filled = true;
                break;
            }
            z_init = ++z;
        }

        #pragma omp parallel for schedule(static, 1)
        for (unsigned int i = 0; i < cell_fifo.size(); ++i) { // parallelize ?
            auto& cell = cell_fifo[i];
            if (!get<2>(cell)->is_cached()) {
                get<2>(cell)->compute();
                if (cache_enabled) {
                    get<2>(cell)->save(cache_context, get<1>(cell));
                }
            }
        }

        for (unsigned int i = 0; i < cell_fifo.size(); ++i) {
            auto& cell = cell_fifo[i];
            get<2>(cell)->extract(index_variant, index_bp, region_start_bp, region_stop_bp, result);
            if (result.last_j >= 0) {
                result.last_cell = get<0>(cell);
                filled = true;
                break;
            }
            if (result.n_correlations >= result.limit) {
                if (i < cell_fifo.size() - 1) {
                    result.last_cell = get<0>(cell_fifo[i + 1]);
                    result.last_j = 0;
                }
                filled = true;
                break;
            }
        }
    }

//    while (z <= z_max) {
//        string key = make_cell_cache_key(cache_key, samples_name, correlation_type, region_chromosome, z);
//        from_morton_code(z, i, j);
//        shared_ptr<Cell> cell = factory.create(correlation_type, i, j);
//        if (cache_enabled) {
//            cell->load(cache_context, key);
//        }
//        cell->segment_i = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_i(), segments);
//        if (!cell->is_diagonal()) {
//            cell->segment_j = load_segment(raw_it->second, store, samples_name, cell->is_cached(), region_chromosome, cell->get_j(), segments);
//        }
//        if (!cell->is_cached()) {
//            cell->compute();
//            if (cache_enabled) {
//                cell->save(cache_context, key);
//            }
//        }
//        cell->extract(index_variant, index_bp, region_start_bp, region_stop_bp, result);
//        if (result.last_j >= 0) {
//            result.last_cell = z;
//            break;
//        }
//        z = get_next_z(segment_index, segment_i, segment_j, z_min, z_max, ++z);
//        if (result.n_correlations >= result.limit) {
//            if (z <= z_max) {
//                result.last_cell = z;
//                result.last_j = 0;
//            }
//            break;
//        }
//    }

    //  Before returning the result, fill out full variant names, chromosomes, and positions
    if (result.raw_index_variant.second >= 0) {
        auto segment = segments[result.raw_index_variant.first];
        result.data.index_variant = segment->get_name(result.raw_index_variant.second);
        result.data.index_chromosome = segment->get_chromosome();
        result.data.index_position = segment->get_position(result.raw_index_variant.second);
    }

    if (result.raw_variants.size() > 0) {
        result.data.variants.reserve(result.raw_variants.size());
        result.data.chromosomes.reserve(result.raw_variants.size());
        result.data.positions.reserve(result.raw_variants.size());
        for (auto&& entry: result.raw_variants) {
            auto segment = segments[entry.first];
            result.data.variants.push_back(segment->get_name(entry.second));
            result.data.chromosomes.push_back(segment->get_chromosome());
            result.data.positions.push_back(segment->get_position(entry.second));
        }
    }

    return true;
}