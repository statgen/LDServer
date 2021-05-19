#ifndef LDSERVER_TYPES_H
#define LDSERVER_TYPES_H

#define  ARMA_DONT_USE_WRAPPER

#include <string>
#include <vector>
#include <algorithm>
#include <regex>
#include <codecvt>
#include <locale>
#include <cmath>
#include <stdexcept>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/writer.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <boost/format.hpp>
#include <armadillo>
#include <msgpack.hpp>
#include <Python.h>
#include <boost/python.hpp>

using namespace std;

class LDServerGenericException : public std::runtime_error {
private:
  std::string secret;
public:
  using std::runtime_error::runtime_error;
  inline LDServerGenericException& set_secret(const std::string& s) { secret = s; return (*this); }
  inline std::string get_secret() const { return secret; }
};

const auto EPACTS_REGEX = regex("(?:chr)?(.+):(\\d+)_?(\\w+)?/?([^_]+)?_?(.*)?");

enum correlation : uint8_t {
    LD_R,
    LD_RSQUARE,
    COV,
    LD_RSQUARE_APPROX
};

/**
 * Enum for column type when reading phenotype files.
 */
enum ColumnType : uint8_t {
  TEXT,
  CATEGORICAL,
  INTEGER,
  FLOAT
};

/**
 * Struct for storing the results of score statistic calculations.
 */
struct ScoreResult {
  string variant;
  double score_stat;
  double pvalue;
  double alt_freq;
  uint64_t position;
  string chrom;

  template<class Archive> void serialize(Archive& ar) { ar(variant, score_stat, pvalue, alt_freq); }
  bool operator==(const ScoreResult& result) const {
      return variant == result.variant;
  }

  bool operator<(const ScoreResult& result) const {
      return position < result.position;
  }
};

/**
 * Common typedefs for shared_ptr
 */
template<class T> using SharedVector = shared_ptr<vector<T>>;
template<typename T> shared_ptr<vector<T>> make_shared_vector(vector<T>& v);
using SharedArmaVec = shared_ptr<arma::vec>;

struct VariantMeta {
    string variant;
    string chromosome;
    string ref;
    string alt;
    uint64_t position;
    string extra;

    /**
     * Construct from an EPACTS-formatted variant. The format looks like:
     *
     *   chrom:pos_ref/alt_extra
     *
     * Where:
     *
     *   chrom - chromosome, can be "chr22" or "22"
     *   pos - position on chromosome
     *   ref - reference allele
     *   alt - alternate allele
     *   extra - usually rsID or sequence identifier, can be any simple string
     *
     * @param variant
     */
    VariantMeta(const string& variant) {
        match_results<string::const_iterator> results;
        regex_search(variant, results, EPACTS_REGEX);

        this->variant = results[0];
        this->chromosome = results[1];
        this->position = stoull(results[2]);
        this->ref = results[3];
        this->alt = results[4];
        this->extra = results[5];
    }

    VariantMeta(const string& variant, const string& chromosome, const string& ref, const string& alt, uint64_t position) : variant(variant), chromosome(chromosome), ref(ref), alt(alt), position(position) {}

    VariantMeta(const string& chromosome, const string& ref, const string& alt, uint64_t position) : chromosome(chromosome), ref(ref), alt(alt), position(position) {
      // These objects store variant name internally as EPACTS format.
      this->variant.append(chromosome)
        .append(":")
        .append(to_string(position))
        .append("_")
        .append(ref)
        .append("/")
        .append(alt);
    }

    bool operator==(VariantMeta const& result) const { // needed by boost.python
        return variant == result.variant;
    }

    bool operator<(VariantMeta const& other) const {
        return position < other.position;
    }

    string as_epacts() const {
      return boost::str(boost::format("%s:%s_%s/%s") % chromosome % position % ref % alt);
    }

    string as_colons() const {
      return boost::str(boost::format("%s:%s:%s:%s") % chromosome % position % ref % alt);
    }
};

struct VariantFrequency {
    string variant;
    string chromosome;
    string ref;
    string alt;
    uint64_t position;
    float ref_af;
    float alt_af;
    VariantFrequency(const string& variant, const string& chromosome, const string& ref, const string& alt, uint64_t position, float ref_af, float alt_af):
            variant(variant), chromosome(chromosome), ref(ref), alt(alt), position(position), ref_af(ref_af), alt_af(alt_af) {}
    bool operator==(VariantFrequency const& result) const { // needed by boost.python
        return variant.compare(result.variant) == 0;
    }
};

/*
 * Struct to store the LD results as upper triangular matrix (for region LD query). This struct is serialized to MessagePack and JSON as a dictionary.
 * All variants are sorted by their position. The order is guaranteed by LDQueryResult.get_variant_index() and LDQueryResult.add_correlation() methods.
 * */
struct LDQueryResultMatrix {
    vector<string> variants;
    vector<string> chromosomes;
    vector<uint64_t> positions;
    vector<int32_t> offsets;
    vector<vector<double>> correlations;
    void clear() {
        variants.clear();
        chromosomes.clear();
        positions.clear();
        offsets.clear();
        correlations.clear();
    }
    MSGPACK_DEFINE_MAP(variants, chromosomes, positions, offsets, correlations) // tell msgpack what fields do we want to put into the JSON-like map
};

/*
 * Struct to store the LD results as a vector (for single variant LD query). This struct is serialized to MessagePack and JSON as a dictionary.
 * All variants are sorted by their position. The order is guaranteed by SingleVariantLDQueryResult.get_variant_index() and SingleVariantLDQueryResult.add_correlation() methods.
 * */
struct LDQueryResultVector {
    string index_variant;
    string index_chromosome;
    uint64_t index_position;
    vector<string> variants;
    vector<string> chromosomes;
    vector<uint64_t> positions;
    vector<double> correlations;
    void clear() {
        index_variant.clear();
        index_chromosome.clear();
        index_position = 0u;
        variants.clear();
        chromosomes.clear();
        positions.clear();
        correlations.clear();
    }
    LDQueryResultVector(): index_variant(""), index_chromosome(""), index_position(0u) {}
    MSGPACK_DEFINE_MAP(index_variant, index_chromosome, index_position, variants, chromosomes, positions, correlations) // tell msgpack what fields do we want to put into the JSON-like map
};

/**
 * Struct to represent LD query result.
 * @param limit This is a constant set in the server config for how many "results" will be returned in one page, in
 *   this case that would be the number of variant pairs.
 * @param last_cell Morton code of cell that should be loaded to begin filling this result object.
 * @param last_i Index of the position to start loading at within the "i" segment of the "last_cell".
 * @param last_j See above, replacing i with j.
 * @param page The current page number. This is incremented by each compute_*_ld() call.
 */
struct LDQueryResult {
    uint32_t limit;
    uint64_t last_cell;
    int last_i;
    int last_j;
    int page;
    unsigned int n_correlations;

    // BEGIN: this is a payload for JSON and MessagePack
    LDQueryResultMatrix data;
    string error;
    string next;
    // END: this is a payload for JSON and MessagePack
    MSGPACK_DEFINE_MAP(data, error, next) // tell msgpack what fields do we want to put into the JSON-like map

    vector<pair<uint64_t, int>> raw_variants; // <segment, index within segment> -- used for temporary collecting all extracted variants

    // defines comparison for <segment, index within segment> pairs to keep variants ordered by chromosomal position
    static bool compare_raw_variants (const pair<uint64_t, int>& v1, const pair<uint64_t, int>& v2) {
        if (v1.first < v2.first) {
            return true;
        }
        return v1.first > v2.first ? false : v1.second < v2.second;
    }

    LDQueryResult(): limit(1000), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0), error(""), next("") { }  // We need this constructor only for MessagePack unit test.
    LDQueryResult(uint32_t page_limit): limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0), error(""), next("") { }

    /**
     * Construct an LDQueryResult
     * @param page_limit
     * @param last This appears to be a string of the form "last_cell:last_i:last_j:page". last_cell is the morton code
     *   of the last cell that was retrieved.
     */
    LDQueryResult(uint32_t page_limit, const string& last) : limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0), error(""), next("") {
        vector<std::string> tokens;
        auto separator = regex(":");
        copy(sregex_token_iterator(last.begin(), last.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
        if (tokens.size() > 3) {
            last_cell = std::stoull(tokens[0]);
            last_i = std::stoi(tokens[1]);
            last_j = std::stoi(tokens[2]);
            page = std::stoi(tokens[3]);
        }
    }

    /**
     * Restrict this result object down to only variants given in the set.
     * @param variants Set of variants.
     */
    template<template <typename...> class C>
    void filter_by_variants(const C<string>& variants) {
        LDQueryResultMatrix new_data;
        unordered_map<uint32_t, uint32_t> index_mapping;
        auto index_mapping_it = index_mapping.end();
        for (unsigned int i = 0u; i < data.variants.size(); ++i) {
            if (variants.find(data.variants[i]) != variants.end()) {
                new_data.variants.push_back(move(data.variants[i]));
                new_data.chromosomes.push_back(move(data.chromosomes[i]));
                new_data.positions.push_back(move(data.positions[i]));
                new_data.offsets.push_back(INT32_MIN);
                new_data.correlations.emplace_back();
                index_mapping.emplace(i, new_data.variants.size() - 1);
            }
        }
        for (unsigned int i1 = 0u; i1 < data.correlations.size(); ++i1) {
            if ((index_mapping_it = index_mapping.find(i1)) == index_mapping.end()) {
                continue;
            }
            auto new_i1 = index_mapping_it->second;
            for (unsigned int j = 0u; j < data.correlations[i1].size(); ++j) {
                auto i2 = j + data.offsets[i1];
                if ((index_mapping_it = index_mapping.find(i2)) == index_mapping.end()) {
                    continue;
                }
                auto new_i2 = index_mapping_it->second;
                new_data.correlations[new_i1].push_back(data.correlations[i1][j]);
                if (new_data.offsets[new_i1] == INT32_MIN) {
                    new_data.offsets[new_i1] = new_i2;
                }
            }
        }
        data = new_data;
    }

    /*
     * If variant is present in "raw_variants", then returns its index there (raw_variants stores variants ordered by chromosomal position).
     * If variant is not present in "raw_variants", then inserts the variant into "raw_variants" and to othe associated data structures (and keeps the ordering).
     * Returns a pair where first element is the index of the variant, and the second element is 0(new variant was not inserted) or 1 (new variant inserted).
     * */
    pair<unsigned int, unsigned int> get_variant(uint64_t segment, int i) {
        auto raw_variant = pair<uint64_t, int>(segment, i);
        unsigned int index = distance(raw_variants.begin(), upper_bound(raw_variants.begin(), raw_variants.end(), raw_variant, compare_raw_variants));
        if ((index > 0) && (raw_variants[index - 1].first == segment) && (raw_variants[index - 1].second == i)) {
            return make_pair(index - 1, 0);
        }
        raw_variants.emplace(raw_variants.begin() + index, move(raw_variant));
        data.offsets.emplace(data.offsets.begin() + index, INT32_MIN);
        data.correlations.emplace(data.correlations.begin() + index, vector<double>());
        for (auto &&o: data.offsets) {
            if ((o != INT32_MIN) && (index <= o)) {
                ++o;
            }
        }
        return make_pair(index, 1);
    }

    // Adds new variant to the result and return its index in sorted data.variants array. If the variant already exists, then does nothing and return its index.
    pair<unsigned int, bool> get_variant_index(const string& variant, const string& chromosome, uint64_t position) {
        unsigned int index = distance(data.positions.begin(), upper_bound(data.positions.begin(), data.positions.end(), position));
        int cmp = 0;
        while ((index > 0) && (data.positions[index - 1] == position)) { // if position in front matches, then we "slide left" comparing variant names.
            cmp = data.variants[index - 1].compare(variant);
            if (cmp == 0) {
                return make_pair(index - 1, false);
            } else if (cmp > 0) {
                --index;
            } else {
                break;
            }
        }
        data.positions.emplace(data.positions.begin() + index, position);
        data.chromosomes.emplace(data.chromosomes.begin() + index, chromosome);
        data.variants.emplace(data.variants.begin() + index, variant);
        data.offsets.emplace(data.offsets.begin() + index, INT32_MIN);
        data.correlations.emplace(data.correlations.begin() + index, vector<double>());
        if (data.offsets.size() - index != 1) {
            for (auto &&o: data.offsets) {
                if ((o != INT32_MIN) && (index <= o)) {
                    ++o;
                }
            }
        }
        return make_pair(index, true);
    }

    tuple<unsigned int, unsigned int, unsigned int> get_variants_range(uint64_t segment, int first_i, int last_i) {
        auto n = last_i - first_i;
        auto last_raw_variant = pair<uint64_t, int>(segment, last_i - 1);
        unsigned int last_index = distance(raw_variants.begin(), upper_bound(raw_variants.begin(), raw_variants.end(), last_raw_variant));
        if (last_index == 0) {
            // all exisiting positions are greater than the last position to insert, so we append everything to the beggining and increment all offsets by number of new variants.
            vector<pair<uint64_t, int>> new_raw_variants(n);
            for (unsigned int i = 0; i < new_raw_variants.size(); ++i) {
                new_raw_variants[i].first = segment;
                new_raw_variants[i].second = first_i + i;
            }
            raw_variants.insert(raw_variants.begin(), new_raw_variants.begin(), new_raw_variants.end());
            data.offsets.insert(data.offsets.begin(), n, INT32_MIN);
            data.correlations.insert(data.correlations.begin(), n, move(vector<double>()));
            for (auto &&o: data.offsets) {
                if (o != INT32_MIN) {
                    o += n;
                }
            }
            return make_tuple(0, 0, n);
        }
        auto first_raw_variant = pair<uint64_t, int>(segment, first_i);
        unsigned int first_index = distance(raw_variants.begin(), lower_bound(raw_variants.begin(), raw_variants.end(), first_raw_variant));
        if (first_index == raw_variants.size()) {
            // all existing positions are smaller than the first position to insert, so we append everything to the end.
            vector<pair<uint64_t, int>> new_raw_variants(n);
            for (unsigned int i = 0; i < new_raw_variants.size(); ++i) {
                new_raw_variants[i].first = segment;
                new_raw_variants[i].second = first_i + i;
            }
            raw_variants.insert(raw_variants.end(), new_raw_variants.begin(), new_raw_variants.end());
            data.offsets.insert(data.offsets.end(), n, INT32_MIN);
            data.correlations.insert(data.correlations.end(), n, move(vector<double>()));
            return make_tuple(raw_variants.size() - n, raw_variants.size() - n, n);
        }
        if (last_index - first_index < n) {
            auto n_new = n - (last_index - first_index);
            vector<pair<uint64_t, int>> new_raw_variants(n_new);
            unsigned int inserted_from = 0u;
            if ((last_index == raw_variants.size()) && ((raw_variants[last_index - 1].first != last_raw_variant.first) || (raw_variants[last_index - 1].second != last_raw_variant.second))) {
                // overlap start
                for (unsigned int i = 0u; i < new_raw_variants.size(); ++i) {
                    new_raw_variants[i].first = segment;
                    new_raw_variants[i].second = first_i + (n - n_new) + i;
                }
                inserted_from = raw_variants.size();
                raw_variants.insert(raw_variants.end(), new_raw_variants.begin(), new_raw_variants.end());
                data.offsets.insert(data.offsets.end(), n_new, INT32_MIN);
                data.correlations.insert(data.correlations.end(), n_new, move(vector<double>()));
            } else {
                // overlap end
                for (unsigned int i = 0u; i < new_raw_variants.size(); ++i) {
                    new_raw_variants[i].first = segment;
                    new_raw_variants[i].second = first_i + i;
                }
                inserted_from = first_index;
                raw_variants.insert(raw_variants.begin() + first_index, new_raw_variants.begin(), new_raw_variants.end());
                data.offsets.insert(data.offsets.begin() + first_index, n_new, INT32_MIN);
                data.correlations.insert(data.correlations.begin() + first_index, n_new, move(vector<double>()));
                for (auto &&o: data.offsets) {
                    if ((o != INT32_MIN) && (first_index <= o)) {
                        o += n_new;
                    }
                }
            }
            return make_tuple(first_index, inserted_from, n_new);
        }
        return make_tuple(first_index, 0, 0);
    }

    // Adds new correlation value between variants at index1 and index2. Must be called only once for a pair i.e. there is no check if the correlation value already exists.
    void add_correlation(uint32_t index1, uint32_t index2, double value) {
        if (index1 > index2) {
            auto temp = index1;
            index1 = index2;
            index2 = temp;
        }
        data.correlations[index1].push_back(value);
        if (data.offsets[index1] == INT32_MIN) {
            data.offsets[index1] = index2;
        }
        ++n_correlations;
    }

    void add_correlations(uint32_t index1, uint32_t start_index2, const float* values, unsigned int n) {
        if (index1 <= start_index2) {
            data.correlations[index1].insert(data.correlations[index1].end(), values, values + n);
            if (data.offsets[index1] == INT32_MIN) {
                data.offsets[index1] = start_index2;
            }
            n_correlations += n;
        } else {
            for (unsigned int i = 0u; i < n; ++i) {
                data.correlations[start_index2 + i].insert(data.correlations[start_index2 + i].end(), values, values + 1);
                if (data.offsets[start_index2 + i] == INT32_MIN) {
                    data.offsets[start_index2 + i] = index1;
                }
                ++values;
                ++n_correlations;
            }
        }
    }

    bool has_next() const {
        return ((last_i >= 0) || (last_j >= 0));
    }
    bool is_last() const {
        return ((page > 0) && (last_i < 0) && (last_j < 0));
    }

    /**
     * Construct a string that represents where to begin loading from. The format is:
     *   "last_cell:last_i:last_j:page",
     * where last_cell is the morton code of the cell to load, and last_i/j are the
     * segment indexes.
     *
     * This function is called by LDQueryResult::get_json() and LDQueryResult::get_messagepack() to provide as part of the "next URL".
     *
     * */
    void set_next(const string& url) {
        if (has_next()) {
            next = url + "&last=" + to_string(last_cell) + ":" + to_string(last_i) + ":" + to_string(last_j) + ":" + to_string(page);
        } else {
            next = "";
        }
    }
    void clear_data() {
        n_correlations = 0u;
        data.clear();
        raw_variants.clear();
        error = "";
        next = "";
    }
    void clear_last() {
        last_cell = 0u;
        last_i = last_j = -1;
    }
    void erase() {
        clear_data();
        clear_last();
        page = 0;
    }
    string get_json(const string& url, int precision = 0) {
        set_next(url);

        rapidjson::StringBuffer strbuf;
        rapidjson::Writer<rapidjson::StringBuffer> writer(strbuf);

        if (precision > 0) {
            writer.SetMaxDecimalPlaces(precision);
        }

//        auto start = std::chrono::system_clock::now();
        writer.StartObject();
        writer.Key("data");
        writer.StartObject();
        writer.Key("variants");
        writer.StartArray();
        for (auto&& v: data.variants) {
            writer.String(v.c_str());
        }
        writer.EndArray();
        writer.Key("chromosomes");
        writer.StartArray();
        for (auto&& v: data.chromosomes) {
            writer.String(v.c_str());
        }
        writer.EndArray();
        writer.Key("positions");
        writer.StartArray();
        for (auto&& v: data.positions) {
            writer.Uint64(v);
        }
        writer.EndArray();
        writer.Key("offsets");
        writer.StartArray();
        for (auto&& v: data.offsets) {
            writer.Int(v);
        }
        writer.EndArray();
        writer.Key("correlations");
        writer.StartArray();
        if (precision > 0) {
            double d = pow(10.0, precision);
            for (auto &&v_correlations : data.correlations) {
                writer.StartArray();
                for (auto &&v: v_correlations) {
                    if (std::isnan(v)) {
                        writer.Null();
                    } else {
                        writer.Double(round(v * d) / d);
                    }
                }
                writer.EndArray();
            }
        } else {
            for (auto &&v_correlations : data.correlations) {
                writer.StartArray();
                for (auto &&v: v_correlations) {
                    if (std::isnan(v)) {
                        writer.Null();
                    } else {
                        writer.Double(v);
                    }
                }
                writer.EndArray();
            }
        }
        writer.EndArray();
        writer.EndObject();
        writer.Key("error");
        writer.String(error.c_str());
        writer.Key("next");
        writer.String(next.c_str());
        writer.EndObject();
//        auto end = std::chrono::system_clock::now();
//        std::chrono::duration<double>  elapsed = end - start;
//        std::cout << "Writing JSON to string elapsed time: " << elapsed.count() << " s\n";
        return strbuf.GetString();
    }

    // Only for unit tests in C/C++
    pair<char*, unsigned int> get_messagepack(const string& url) {
        set_next(url);
//        auto start = std::chrono::system_clock::now();
        msgpack::sbuffer msgpack_buffer;
        msgpack::pack(msgpack_buffer, *this);
//        auto end = std::chrono::system_clock::now();
//        std::chrono::duration<double>  elapsed = end - start;
//        std::cout << "Packing msgpack elapsed time: " << elapsed.count() << " s\n";
        unsigned int size = msgpack_buffer.size();
        return make_pair(msgpack_buffer.release(), size);
    }

    PyObject* get_messagepack_py(const string& url) {
        auto msgpack = get_messagepack(url);
        PyObject* msgpack_py = PyBytes_FromStringAndSize(msgpack.first, msgpack.second);
        free(msgpack.first);
//        boost::python::object memoryView(boost::python::handle<>(PyMemoryView_FromMemory(msgpack_buffer.data(), msgpack_buffer.size(), PyBUF_READ)));
//        boost::python::object memoryView(boost::python::handle<>(PyBytes_FromStringAndSize(msgpack_buffer.data(), msgpack_buffer.size())));
        return msgpack_py;
    }
};


struct SingleVariantLDQueryResult {
    uint32_t limit;
    uint64_t last_cell;
    int last_i;
    int last_j;
    int page;
    unsigned int n_correlations;

    // BEGIN: this is a payload for JSON and MessagePack
    LDQueryResultVector data;
    string error;
    string next;
    // END: this is a payload for JSON and MessagePack
    MSGPACK_DEFINE_MAP(data, error, next) // tell msgpack what fields do we want to put into the JSON-like map

    pair<uint64_t, int> raw_index_variant; // <segment, index within segment> -- used for temporary storing index variant
    vector<pair<uint64_t, int>> raw_variants; // <segment, index within segment> -- used for temporary collecting all extracted variants

    // defines comparison for <segment, index within segment> pairs to keep variants ordered by chromosomal position
    static bool compare_raw_variants (const pair<uint64_t, int>& v1, const pair<uint64_t, int>& v2) {
        if (v1.first < v2.first) {
            return true;
        }
        return v1.first > v2.first ? false : v1.second < v2.second;
    }

    SingleVariantLDQueryResult(): limit(1000), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0), error(""), next(""), raw_index_variant(0, -1) { }  // We need this constructor only for MessagePack unit test.
    SingleVariantLDQueryResult(uint32_t page_limit): limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0), error(""), next(""), raw_index_variant(0, -1) { }

    /**
     * Construct an SingleVariantLDQueryResult
     * @param page_limit
     * @param last This appears to be a string of the form "last_cell:last_i:last_j:page". last_cell is the morton code
     *   of the last cell that was retrieved.
     */
    SingleVariantLDQueryResult(uint32_t page_limit, const string& last) : limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0), error(""), next(""), raw_index_variant(0, -1) {
        vector<std::string> tokens;
        auto separator = regex(":");
        copy(sregex_token_iterator(last.begin(), last.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
        if (tokens.size() > 3) {
            last_cell = std::stoull(tokens[0]);
            last_i = std::stoi(tokens[1]);
            last_j = std::stoi(tokens[2]);
            page = std::stoi(tokens[3]);
        }
    }

    // returns 1 - if new variant was set; returns 0 if no new index variant was set
    int get_index_variant(uint64_t segment, int i) {
        if (raw_index_variant.second > 0) {
            return 0;
        }
        raw_index_variant.first = segment;
        raw_index_variant.second = i;
        return i;
    }

//    /*
//     * If variant is present in "raw_variants", then returns its index there (raw_variants stores variants ordered by chromosomal position).
//     * If variant is not present in "raw_variants", then inserts the variant into "raw_variants" and to othe associated data structures (and keeps the ordering).
//     * Returns a pair where first element is the index of the variant, and the second element is 0(new variant was not inserted) or 1 (new variant inserted).
//     * */
//    pair<unsigned int, unsigned int> get_variant(uint64_t segment, int i) {
//        auto raw_variant = pair<uint64_t, int>(segment, i);
//        unsigned int index = distance(raw_variants.begin(), upper_bound(raw_variants.begin(), raw_variants.end(), raw_variant, compare_raw_variants));
//        if ((index > 0) && (raw_variants[index - 1].first == segment) && (raw_variants[index - 1].second == i)) {
//            return make_pair(index - 1, 0);
//        }
//        raw_variants.emplace(raw_variants.begin() + index, move(raw_variant));
////        data.offsets.emplace(data.offsets.begin() + index, INT32_MIN);
////        data.correlations.emplace(data.correlations.begin() + index, vector<double>());
////        for (auto &&o: data.offsets) {
////            if ((o != INT32_MIN) && (index <= o)) {
////                ++o;
////            }
////        }
//        return make_pair(index, 1);
//    }

    tuple<unsigned int, unsigned int, unsigned int> get_variants_range(uint64_t segment, int first_i, int last_i) {
        auto n = last_i - first_i;
        auto last_raw_variant = pair<uint64_t, int>(segment, last_i - 1);
        unsigned int last_index = distance(raw_variants.begin(), upper_bound(raw_variants.begin(), raw_variants.end(), last_raw_variant));
        if (last_index == 0) {
            // all exisiting positions are greater than the last position to insert, so we append everything to the beggining and increment all offsets by number of new variants.
            vector<pair<uint64_t, int>> new_raw_variants(n);
            for (unsigned int i = 0; i < new_raw_variants.size(); ++i) {
                new_raw_variants[i].first = segment;
                new_raw_variants[i].second = first_i + i;
            }
            raw_variants.insert(raw_variants.begin(), new_raw_variants.begin(), new_raw_variants.end());
            return make_tuple(0, 0, n);
        }
        auto first_raw_variant = pair<uint64_t, int>(segment, first_i);
        unsigned int first_index = distance(raw_variants.begin(), lower_bound(raw_variants.begin(), raw_variants.end(), first_raw_variant));
        if (first_index == raw_variants.size()) {
            // all existing positions are smaller than the first position to insert, so we append everything to the end.
            vector<pair<uint64_t, int>> new_raw_variants(n);
            for (unsigned int i = 0; i < new_raw_variants.size(); ++i) {
                new_raw_variants[i].first = segment;
                new_raw_variants[i].second = first_i + i;
            }
            raw_variants.insert(raw_variants.end(), new_raw_variants.begin(), new_raw_variants.end());
            return make_tuple(raw_variants.size() - n, raw_variants.size() - n, n);
        }
        if (last_index - first_index < n) {
            auto n_new = n - (last_index - first_index);
            vector<pair<uint64_t, int>> new_raw_variants(n_new);
            unsigned int inserted_from = 0u;
            if ((last_index == raw_variants.size()) && ((raw_variants[last_index - 1].first != last_raw_variant.first) || (raw_variants[last_index - 1].second != last_raw_variant.second))) {
                // overlap start
                for (unsigned int i = 0u; i < new_raw_variants.size(); ++i) {
                    new_raw_variants[i].first = segment;
                    new_raw_variants[i].second = first_i + (n - n_new) + i;
                }
                inserted_from = raw_variants.size();
                raw_variants.insert(raw_variants.end(), new_raw_variants.begin(), new_raw_variants.end());
            } else {
                // overlap end
                for (unsigned int i = 0u; i < new_raw_variants.size(); ++i) {
                    new_raw_variants[i].first = segment;
                    new_raw_variants[i].second = first_i + i;
                }
                inserted_from = first_index;
                raw_variants.insert(raw_variants.begin() + first_index, new_raw_variants.begin(), new_raw_variants.end());
            }
            return make_tuple(first_index, inserted_from, n_new);
        }
        return make_tuple(first_index, 0, 0);
    }

//    // Adds new correlation value between variants at index1 and index2. Must be called only once for a pair i.e. there is no check if the correlation value already exists.
//    void add_correlation(uint32_t index1, uint32_t index2, double value) {
////        if (index1 > index2) {
////            auto temp = index1;
////            index1 = index2;
////            index2 = temp;
////        }
////        data.correlations[index1].push_back(value);
////        if (data.offsets[index1] == INT32_MIN) {
////            data.offsets[index1] = index2;
////        }
//        ++n_correlations;
//    }

    void add_correlations(uint32_t start_index, const float* values, unsigned int n) {
        data.correlations.insert(data.correlations.begin() + start_index, values, values + n);
        n_correlations += n;
    }

    bool has_next() const {
        return ((last_i >= 0) || (last_j >= 0));
    }

    bool is_last() const {
        return ((page > 0) && (last_i < 0) && (last_j < 0));
    }

    /**
     * Construct a string that represents where to begin loading from. The format is:
     *   "last_cell:last_i:last_j:page",
     * where last_cell is the morton code of the cell to load, and last_i/j are the
     * segment indexes.
     *
     * This function is called by LDQueryResult::get_json() and LDQueryResult::get_messagepack() to provide as part of the "next URL".
     *
     * */
    void set_next(const string& url) {
        if (has_next()) {
            next = url + "&last=" + to_string(last_cell) + ":" + to_string(last_i) + ":" + to_string(last_j) + ":" + to_string(page);
        } else {
            next = "";
        }
    }

    void clear_data() {
        n_correlations = 0u;
        data.clear();
        raw_index_variant.first = 0u;
        raw_index_variant.second = -1;
        raw_variants.clear();
        error = "";
        next = "";
    }

    void clear_last() {
        last_cell = 0u;
        last_i = last_j = -1;
    }

    void erase() {
        clear_data();
        clear_last();
        page = 0;
    }

    string get_json(const string& url, int precision = 0) {
        set_next(url);

        rapidjson::StringBuffer strbuf;
        rapidjson::Writer<rapidjson::StringBuffer> writer(strbuf);

        if (precision > 0) {
            writer.SetMaxDecimalPlaces(precision);
        }

//        auto start = std::chrono::system_clock::now();
        writer.StartObject();
        writer.Key("data");
        writer.StartObject();
        writer.Key("index_variant");
        if (!data.index_variant.empty()) {
            writer.String(data.index_variant.c_str());
        } else {
            writer.Null();
        }
        writer.Key("index_chromosome");
        if (!data.index_chromosome.empty()) {
            writer.String(data.index_chromosome.c_str());
        } else {
            writer.Null();
        }
        writer.Key("index_position");
        if (data.index_position > 0) {
            writer.Uint64(data.index_position);
        } else {
            writer.Null();
        }
        writer.Key("variants");
        writer.StartArray();
        for (auto&& v: data.variants) {
            writer.String(v.c_str());
        }
        writer.EndArray();
        writer.Key("chromosomes");
        writer.StartArray();
        for (auto&& v: data.chromosomes) {
            writer.String(v.c_str());
        }
        writer.EndArray();
        writer.Key("positions");
        writer.StartArray();
        for (auto&& v: data.positions) {
            writer.Uint64(v);
        }
        writer.EndArray();
        writer.Key("correlations");
        writer.StartArray();
        if (precision > 0) {
            double d = pow(10.0, precision);
            for (auto& v : data.correlations) {
                if (std::isnan(v)) {
                    writer.Null();
                } else {
                    writer.Double(round(v * d) / d);
                }
            }
        } else {
            for (auto& v : data.correlations) {
                if (std::isnan(v)) {
                    writer.Null();
                } else {
                    writer.Double(v);
                }
            }
        }
        writer.EndArray();
        writer.EndObject();
        writer.Key("error");
        writer.String(error.c_str());
        writer.Key("next");
        writer.String(next.c_str());
        writer.EndObject();
//        auto end = std::chrono::system_clock::now();
//        std::chrono::duration<double>  elapsed = end - start;
//        std::cout << "Writing JSON to string elapsed time: " << elapsed.count() << " s\n";
        return strbuf.GetString();
    }

    // Only for unit tests in C/C++
    pair<char*, unsigned int> get_messagepack(const string& url) {
        set_next(url);
//        auto start = std::chrono::system_clock::now();
        msgpack::sbuffer msgpack_buffer;
        msgpack::pack(msgpack_buffer, *this);
//        auto end = std::chrono::system_clock::now();
//        std::chrono::duration<double>  elapsed = end - start;
//        std::cout << "Packing msgpack elapsed time: " << elapsed.count() << " s\n";
        unsigned int size = msgpack_buffer.size();
        return make_pair(msgpack_buffer.release(), size);
    }

    PyObject* get_messagepack_py(const string& url) {
        auto msgpack = get_messagepack(url);
        PyObject* msgpack_py = PyBytes_FromStringAndSize(msgpack.first, msgpack.second);
        free(msgpack.first);
//        boost::python::object memoryView(boost::python::handle<>(PyMemoryView_FromMemory(msgpack_buffer.data(), msgpack_buffer.size(), PyBUF_READ)));
//        boost::python::object memoryView(boost::python::handle<>(PyBytes_FromStringAndSize(msgpack_buffer.data(), msgpack_buffer.size())));
        return msgpack_py;
    }
};

/**
 * Struct to represent results of a score stat query.
 * @param limit The maximum number of results to be returned.
 * @param last_i The last index returned from within a segment.
 * @param last_seg The last segment that was in use
 * @param page Which page are we on
 * @param data Vector of results.
 */
struct ScoreStatQueryResult {
  uint64_t limit;
  int64_t last_i;
  int64_t last_seg;
  uint64_t page;
  vector<ScoreResult> data;
  double sigma2; // assume all score stats calculated against the same phenotype
  double nsamples;

  ScoreStatQueryResult(uint32_t page_limit): limit(page_limit), last_i(-1), last_seg(0), page(0) {
      data.reserve(page_limit);
  }

  ScoreStatQueryResult(uint32_t page_limit, const string& last) : limit(page_limit), last_i(-1), last_seg(0), page(0) {
      data.reserve(page_limit);
      vector<std::string> tokens;
      auto separator = regex(":");
      copy(sregex_token_iterator(last.begin(), last.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
      if (tokens.size() > 3) {
          last_seg = std::stoll(tokens[0]);
          last_i = std::stoll(tokens[1]);
          page = std::stoull(tokens[2]);
      }
  }

  void sort_by_variant() {
      std::sort(data.begin(), data.end(), [](const ScoreResult& a, const ScoreResult& b) {
          return a.position < b.position;
      });
  }

  /**
   * Restrict this result object down to only variants given in the set.
   * @param variants Set of variants.
   */
  template<template <typename...> class C>
  void filter_by_variants(const C<string>& variants) {
    vector<ScoreResult> new_data;
    copy_if(data.begin(), data.end(), back_inserter(new_data), [&](const ScoreResult& p) -> bool {
      return variants.find(p.variant) != variants.end();
    });
    data = new_data;
  }

  bool has_next() const {
      return last_i >= 0;
  }

  bool is_last() const {
      return ((page > 0) && (last_i < 0));
  }

  string get_last() const {
      if (has_next()) {
          return string(to_string(last_seg) + ":" + to_string(last_i) + ":" + to_string(page));
      }
      return string("");
  }

  void clear_data() {
      data.clear();
      sigma2 = numeric_limits<double>::quiet_NaN();
      nsamples = numeric_limits<double>::quiet_NaN();
  }

  void clear_last() {
      last_seg = 0u;
      last_i = -1;
  }

  void erase() {
      clear_data();
      clear_last();
      sigma2 = numeric_limits<double>::quiet_NaN();
      nsamples = numeric_limits<double>::quiet_NaN();
      page = 0;
  }

//  struct ScoreResult {
//    string variant;
//    double score_stat;
//    double sigma2;
//    double pvalue;
//    double alt_freq;
//  };
  string get_json(const string& url) {
      rapidjson::Document document;
      document.SetObject();
      rapidjson::Document::AllocatorType& allocator = document.GetAllocator();

      rapidjson::Value data(rapidjson::kObjectType);

      rapidjson::Value variant(rapidjson::kArrayType);
      rapidjson::Value alt_freq(rapidjson::kArrayType);
      rapidjson::Value pvalue(rapidjson::kArrayType);
      rapidjson::Value score_stat(rapidjson::kArrayType);

      for (auto&& p: this->data) {
          variant.PushBack(rapidjson::StringRef(p.variant.c_str()), allocator);
//          score_stat.PushBack(p.score_stat, allocator);
//          sigma2.PushBack(p.sigma2, allocator);

          if (std::isnan(p.pvalue)) { // nan is not allowed by JSON, so we replace it with null
              pvalue.PushBack(rapidjson::Value(), allocator);
          }
          else {
              pvalue.PushBack(p.pvalue, allocator);
          }

          if (std::isnan(p.score_stat)) { // nan is not allowed by JSON, so we replace it with null
              score_stat.PushBack(rapidjson::Value(), allocator);
          }
          else {
              score_stat.PushBack(p.score_stat, allocator);
          }
      }

      data.AddMember("variant", variant, allocator);
      data.AddMember("alt_freq", alt_freq, allocator);
      data.AddMember("pvalue", pvalue, allocator);
      data.AddMember("score_stat", score_stat, allocator);

      if (std::isnan(sigma2)) { // nan is not allowed by JSON, so we replace it with null
        data.AddMember("sigma2", rapidjson::Value(), allocator);
      }
      else {
        data.AddMember("sigma2", sigma2, allocator);
      }

      if (std::isnan(nsamples)) { // nan is not allowed by JSON, so we replace it with null
        data.AddMember("n_samples", rapidjson::Value(), allocator);
      }
      else {
        data.AddMember("n_samples", nsamples, allocator);
      }

      document.AddMember("data", data, allocator);
      document.AddMember("error", rapidjson::Value(), allocator);

      if (is_last()) {
          document.AddMember("next", rapidjson::Value(), allocator);
      }
      else {
          rapidjson::Value next;
          string link = url + "&last=" + get_last();
          next.SetString(link.c_str(), allocator); // by providing allocator we make a copy
          document.AddMember("next", next, allocator);
      }

      rapidjson::StringBuffer strbuf;
      rapidjson::Writer<rapidjson::StringBuffer> writer(strbuf);
      if (!document.Accept(writer)) {
          throw runtime_error("Error while saving to JSON");
      }

      return strbuf.GetString();
  }
};

#endif
