#ifndef LDSERVER_LDSERVER_H
#define LDSERVER_LDSERVER_H

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
#include "Raw.h"
#include "Morton.h"
#include "Segment.h"
#include "Cell.h"
#include "Types.h"

using namespace std;

class LDServer {
private:
    unordered_map<string, vector<string>> samples;
    unordered_map<string, shared_ptr<Raw>> raw;

    bool cache_enabled;
    string cache_hostname;
    int cache_port;
    redisContext* cache_context;

    static void parse_variant(const string& variant, string& chromosome, uint64_t& position, string& ref_allele, string& alt_allele);

public:
    static const string ALL_SAMPLES_KEY;

    LDServer();
    virtual ~LDServer();

    void set_file(const string& file);
    void set_samples(const string& name, const vector<string>& samples);
    void enable_cache(const string& hostname, int port);
    void disable_cache();

    bool compute_region_ld(const string& region_chromosome, uint64_t region_start_bp, uint64_t region_stop_bp, struct LDQueryResult& result, const string& samples_name = ALL_SAMPLES_KEY) const;
    bool compute_variant_ld(const string& index_variant, const string& region_chromosome, uint64_t region_start_bp, uint64_t region_stop_bp, struct LDQueryResult& result, const string& samples_name = ALL_SAMPLES_KEY) const;
};

#endif