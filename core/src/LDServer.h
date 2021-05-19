#ifndef LDSERVER_LDSERVER_H
#define LDSERVER_LDSERVER_H

#define  ARMA_DONT_USE_WRAPPER

#include <string>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <unordered_map>
#include <savvy/reader.hpp>
#include <savvy/armadillo_vector.hpp>
#include <hiredis/hiredis.h>
#include <armadillo>
#include <cmath>
#include <regex>
#include <chrono>
#include "Raw.h"
#include "Morton.h"
#include "Segment.h"
#include "Cell.h"
#include "Types.h"
#include <omp.h>

using namespace std;

class LDServer {
private:
    unordered_map<string, vector<string>> samples;
    unordered_map<string, shared_ptr<Raw>> raw;

    uint32_t segment_size;
    set<uint64_t> allowed_segments;

    bool cache_enabled;
    uint32_t cache_key;
    string cache_hostname;
    int cache_port;
    redisContext* cache_context;

    static void parse_variant(const string& variant, string& chromosome, uint64_t& position, string& ref_allele, string& alt_allele);
    shared_ptr<Segment> load_segment(const shared_ptr<Raw>& raw, genotypes_store store, const string& samples_name, bool only_variants, const std::string& chromosome, uint64_t i, std::map<std::uint64_t, shared_ptr<Segment>>& segments) const;

public:
    static const string ALL_SAMPLES_KEY;

    LDServer(uint32_t segment_size = 1000);
    virtual ~LDServer();

    /**
     * Specify variant positions directly to the server, so that only segments containing these positions
     * are loaded (this is useful for sparse queries with a large range, but few variants.)
     * @param pos Starting position of the variant
     */
    void add_overlap_position(const uint64_t& pos);
    void set_file(const string& file);
    void set_samples(const string& name, const vector<string>& samples);
    void force_samples(const std::string &name, const std::vector<std::string> &samples);
    void enable_cache(uint32_t cache_key, const string& hostname, int port);
    void disable_cache();

    static string make_cell_cache_key(uint32_t cache_key, const string& samples_name, correlation correlation_type, const string& chromosome, uint64_t morton_code);
    static string make_segment_cache_key(uint32_t cache_key, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp);

    vector<string> get_chromosomes();
    uint32_t get_segment_size() const;

    /**
     * Compute LD between all variants in a region.
     * @param region_chromosome
     * @param region_start_bp
     * @param region_stop_bp
     * @param correlation_type
     * @param result
     * @param samples_name
     * @param diagonal Should we compute the diagonal elements? (variance of each variant)
     * @param segments_out Shared vector of segments that can be passed on to the ScoreServer for extra computations.
     * @return
     */
    bool compute_region_ld(const string& region_chromosome, uint64_t region_start_bp, uint64_t region_stop_bp, correlation correlation_type, struct LDQueryResult& result, const string& samples_name = ALL_SAMPLES_KEY, bool diagonal = false, SharedSegmentVector segments_out = nullptr) const;
    bool compute_variant_ld(const string& index_variant, const string& region_chromosome, uint64_t region_start_bp, uint64_t region_stop_bp, correlation correlation_type, struct SingleVariantLDQueryResult& result, const string& samples_name = ALL_SAMPLES_KEY) const;
};

#endif