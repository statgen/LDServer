#include "ScoreServer.h"
#include <stdexcept>

const std::string ScoreServer::ALL_SAMPLES_KEY("ALL");

ScoreServer::ScoreServer(uint32_t segment_size) : segment_size(segment_size), cache_enabled(false), cache_hostname(""), cache_port(0), cache_context(nullptr) {}

ScoreServer::~ScoreServer() {
    if (cache_context != nullptr) {
        redisFree(cache_context);
        cache_context = nullptr;
    }
}

void ScoreServer::set_genotypes_file(const std::string &file, const uint32_t& genotype_dataset_id) {
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
    this->genotype_dataset_id = genotype_dataset_id;
}

void ScoreServer::load_phenotypes_file(const string &path, const ColumnTypeMap &types, size_t nrows,
                                       const string& delim, const string& sample_column, const uint32_t &phenotype_dataset_id) {
    // Release previous phenotypes object.
    phenotypes.reset();

    // Create a new one.
    phenotypes = make_shared<Phenotypes>();

    // Load phenotypes into memory.
    phenotypes->load_file(path, types, nrows, delim, sample_column);

    // Store ID for use in cache key.
    this->phenotype_dataset_id = phenotype_dataset_id;

    #ifndef NDEBUG
    phenotypes->pprint();
    #endif
}

void ScoreServer::set_phenotype(const string& p) {
    // Set the phenotype.
    phenotype = p;
}


shared_ptr<vector<string>> ScoreServer::get_complete_samples(const string& phenotype) const {
  return phenotypes->get_complete_samples(phenotype);
}

void ScoreServer::set_samples(const std::string &name, const std::vector<std::string> &samples) {
    auto e = this->samples.emplace(name, std::vector<std::string>()).first;
    for (auto&& sample: samples) {
        e->second.push_back(sample);
    }
}

void ScoreServer::force_samples(const std::string &name, const std::vector<std::string> &samples) {
    this->samples[name] = samples;
}

/**
 * Enable the cache
 * @param hostname
 * @param port
 */
void ScoreServer::enable_cache(const string& hostname, int port) {
    if (cache_context == nullptr) {
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

void ScoreServer::disable_cache() {
    if (cache_enabled) {
        if (cache_context != nullptr) {
            redisFree(cache_context);
            cache_context = nullptr;
        }
        cache_enabled = false;
    }
}

string ScoreServer::make_segment_cache_key(uint32_t genotype_dataset_id, uint32_t phenotype_dataset_id, const string& phenotype_name, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp) {
    stringstream os(ios::binary | ios::out);
    os.write(reinterpret_cast<const char*>(&genotype_dataset_id), sizeof(genotype_dataset_id));
    os.write(reinterpret_cast<const char*>(&phenotype_dataset_id), sizeof(phenotype_dataset_id));
    os.write(phenotype_name.c_str(), phenotype_name.size());
    os.write(samples_name.c_str(), samples_name.size());
    os.write(chromosome.c_str(), chromosome.size());
    os.write(reinterpret_cast<const char*>(&start_bp), sizeof(start_bp));
    os.write(reinterpret_cast<const char*>(&stop_bp), sizeof(stop_bp));
    os.flush();
    return os.str();
}

void ScoreServer::parse_variant(const std::string& variant, std::string& chromosome, uint64_t& position, std::string& ref_allele, std::string& alt_allele) {
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

vector<string> ScoreServer::get_chromosomes() {
    vector<string> chromosomes;
    for (const auto& x : this->raw) {
        chromosomes.emplace_back(x.first);
    }
    return chromosomes;
}

uint32_t ScoreServer::get_segment_size() const {
    return this->segment_size;
}

/**
 * Compute score statistics on segments passed from the LDServer.
 * @param result
 * @param samples_name
 * @param segments
 * @return
 */
bool ScoreServer::compute_scores(const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct ScoreStatQueryResult& result, const std::string& samples_name, SharedSegmentVector segments) const {
    // Prepare result object
    if (result.is_last()) {
        return false;
    }
    result.clear_data();
    result.page += 1;

    // The LDServer shouldn't pass us an empty vector of segments
    if (segments->empty()) {
        throw std::invalid_argument("Segment vector is empty");
    }

    // Check segments
    // Assumptions:
    //   * All segments are on the same chromosome
    //   * Segments should have had names already loaded by the LDServer
    //const string& region_chromosome = (*segments)[0]->get_chromosome();
    for (auto&& segment : *segments) {
        if (segment->get_chromosome() != region_chromosome) {
            throw std::logic_error("All segments must be of the same chromosome");
        }
        if (!segment->has_names()) {
            throw std::invalid_argument("Error: variant IDs have not been loaded on passed segment to score server");
        }
    }

    // Get the list of samples to use
    auto samples_it = this->samples.find(samples_name);
    if (samples_it == this->samples.end()) { // no such samples - return empty result
        result.clear_last();
        return false;
    }

    // Tell the phenotypes to re-order according to these samples
    phenotypes->reorder(samples_it->second);

    // Find the proper VCF or SAV file to load from
    auto raw_it = raw.find(region_chromosome);
    if (raw_it == raw.end()) { // no such chromosome - return empty result
        result.clear_last();
        return false;
    }

    // This tells the "Raw" object to open a handle to the proper file / chromosome.
    // The final "true" just means retrieve genotype values coded 0/1/2 (additive), which is what is needed
    // when doing covariance/scores.
    raw_it->second->open(region_chromosome, samples_it->second, true);

    // Calculate phenotypic variance
    result.sigma2 = phenotypes->compute_sigma2(this->phenotype);
    result.nsamples = phenotypes->get_nsamples(this->phenotype);

    // Process each segment and run calculations if necessary
    for (int seg_i = 0; seg_i < segments->size(); seg_i++) {
        //TODO: add cache retrieval/storage
        //string segment_key = make_segment_cache_key(genotype_dataset_id, phenotype_dataset_id, phenotype, samples_name, segment->get_chromosome(), segment->get_start_bp(), segment->get_stop_bp());

        // Create a ScoreSegment object from a Segment, moving elements into it
        shared_ptr<ScoreSegment> score_seg = make_shared<ScoreSegment>(std::move(*((*segments)[seg_i])));

        // Load the genotypes if they're not loaded already
        if (!score_seg->has_genotypes()) {
            raw_it->second->load_genotypes(*(score_seg));
        }

        // Perform score statistic calculations
        score_seg->compute_scores(*phenotypes->as_vec(this->phenotype));

        // Extract into results object.
        score_seg->extract(region_start_bp, region_stop_bp, result);

        if (result.last_i >= 0) {
            // Either the segment was only partially iterated over (due to result object limit being reached), or the
            // result object filled up and we needed to stop.
            result.last_seg = seg_i;
            break;
        }

        // If we made it here, segment is done loading.
        // The result object may be full, however.
        if (result.data.size() >= result.limit) {
            if (seg_i+1 < segments->size()) {
                result.last_seg = seg_i + 1;
            }
            break;
        }
    }

    return true;
}