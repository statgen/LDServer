#ifndef LDSERVER_SCORESERVER_H
#define LDSERVER_SCORESERVER_H

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
#include "Segment.h"
#include "ScoreSegment.h"
#include "Types.h"
#include "Phenotypes.h"

using namespace std;

class ScoreServer {
private:
    unordered_map<string, vector<string>> samples;
    unordered_map<string, shared_ptr<Raw>> raw;
    shared_ptr<Phenotypes> phenotypes;
    uint32_t genotype_dataset_id;
    uint32_t phenotype_dataset_id;
    string phenotype;

    uint32_t segment_size;

    bool cache_enabled;
    string cache_hostname;
    int cache_port;
    redisContext* cache_context;

    static void parse_variant(const string& variant, string& chromosome, uint64_t& position, string& ref_allele, string& alt_allele);

public:
    static const string ALL_SAMPLES_KEY;

    explicit ScoreServer(uint32_t segment_size = 1000);
    virtual ~ScoreServer();

    /**
     * Set a file from which to read genotypes. This can be a VCF or Savvy formatted file.
     * @param file
     * @param genotype_dataset_id - integer representing the genotype dataset (IDs stored in sqlite database).
     */
    void set_genotypes_file(const string& file, const uint32_t& genotype_dataset_id);

    /**
     * Set the list of samples to be used. This will set it both when operating on genotypes, and phenotypes.
     * @param name
     * @param samples
     */
    void set_samples(const string& name, const vector<string>& samples);

    void force_samples(const std::string &name, const std::vector<std::string> &samples);

    /**
     * Load phenotypes from a file. This triggers a load of phenotypes into memory. Phenotype files are assumed to
     * typically be quite small. When UKBB gets metabolites, then maybe we will need to start loading only needed
     * phenotypes in a single pass, or switch to alternative storage.
     *
     * These are named "load_" because they load the entire file. The "set_genotypes_file" method above does not,
     * rather it can selectively load genotypes using an on-disk index.
     *
     * @param pedpath
     * @param phenotype_dataset_id (see genotype_dataset_id above)
     */
    void load_phenotypes_file(const string &path, const ColumnTypeMap &types, size_t nrows, const string& delim, const string& sample_column, const uint32_t &phenotype_dataset_id);

    /**
     * Set the phenotype that will be used for score/p-value calculations.
     * @param p
     */
    void set_phenotype(const string& p);

    /**
     * Get list of samples for which the phenotype has complete data.
     */
     shared_ptr<vector<string>> get_complete_samples(const string& phenotype) const;

    /**
     * Functions to enable/disable the cache.
     * Note: this is slightly different than LDServer, where a cache key is not used (these functions only open/close
     * the redis context.)
     * @param hostname
     * @param port
     */
    void enable_cache(const string& hostname, int port);
    void disable_cache();

    /**
     * Create a key to cache a segment of score statistics.
     * @param genotype_dataset_id - ID of genotypes dataset
     * @param phenotype_dataset_id - ID of phenotypes dataset
     * @param phenotype_name - Name of the phenotype
     * @param samples_name - Name of sample subset in use, or "ALL"
     * @param chromosome
     * @param start_bp
     * @param stop_bp
     * @return
     */
    static string make_segment_cache_key(uint32_t genotype_dataset_id, uint32_t phenotype_dataset_id, const string& phenotype_name, const string& samples_name, const string& chromosome, uint64_t start_bp, uint64_t stop_bp);

    vector<string> get_chromosomes();
    uint32_t get_segment_size() const;

    /**
     * Compute score statistics on segments that have already been used by the LDServer.
     * If the segments have cached score statistics, we use them directly.
     * Otherwise, if the segments at least have loaded genotypes, we can calculate the score statistics from them.
     * Finally, if neither is true, we will have to load the genotypes from disk and compute.
     *
     * If the cache is enabled, score statistics and relevant other stats will be stored to the redis cache.
     *
     * @param result
     * @param samples_name
     * @param segments
     * @return
     */
    bool compute_scores(const string& region_chromosome, uint64_t region_start_bp, uint64_t region_stop_bp, struct ScoreStatQueryResult& result, const string& samples_name = ALL_SAMPLES_KEY, SharedSegmentVector segments = nullptr) const;
};

#endif