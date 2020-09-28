#ifndef LDSERVER_TYPES_H
#define LDSERVER_TYPES_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <regex>
#include <codecvt>
#include <locale>
#include <cmath>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/writer.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <armadillo>

using namespace std;

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

    bool operator==(VariantMeta const& result) const { // needed by boost.python
        return variant == result.variant;
    }

    bool operator<(VariantMeta const& other) const {
        return position < other.position;
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

struct VariantsPair {
    string variant1;
    string variant2;
    string chromosome1;
    string chromosome2;
    uint64_t position1;
    uint64_t position2;
    double value;
    VariantsPair() : variant1(""), variant2(""), chromosome1(""), chromosome2(""), position1(0ul), position2(0ul), value(0.0) {}
    VariantsPair(const string& variant1, const string& chromosome1, uint64_t position1, const string& variant2, const string& chromosome2, uint64_t position2, double value):
            variant1(variant1), variant2(variant2), chromosome1(chromosome1), chromosome2(chromosome2), position1(position1), position2(position2), value(value) {}
    bool operator==(VariantsPair const& pair) const { // needed by boost.python
        return (variant1.compare(pair.variant1) == 0 && variant2.compare(pair.variant2) == 0);
    }
};

struct Variant {
    string name;
    string chromosome;
    uint64_t position;
    Variant() : name(""), chromosome(""), position(0u) {}
    Variant(const string& name, const string& chromosome, uint64_t position): name(name), chromosome(chromosome), position(position) {}
    bool operator==(Variant const& other) const { // needed by boost.python
        return (name.compare(other.name) == 0);
    }
    bool operator<(Variant const& other) const {
        if (position == other.position) {
            return name.compare(other.name) < 0;
        }
        return position < other.position;
    }
};

struct Correlation {
    uint32_t variant_idx;
    double value;
    Correlation() : variant_idx(0u), value(0.0) {}
    Correlation(uint32_t variant_idx, double value): variant_idx(variant_idx), value(value) {}
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

    vector<Variant> variants;
    unordered_map<string, uint32_t> index;
    unordered_map<uint32_t, vector<Correlation>> correlations;
    unsigned int n_correlations;

//    vector<VariantsPair> data; // TODO: DT. remove it later


    LDQueryResult(uint32_t page_limit): limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0) {
//        data.reserve(page_limit);
    }

    /**
     * Construct an LDQueryResult
     * @param page_limit
     * @param last This appears to be a string of the form "last_cell:last_i:last_j:page". last_cell is the morton code
     *   of the last cell that was retrieved.
     */
    LDQueryResult(uint32_t page_limit, const string& last) : limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0), n_correlations(0) {
//        data.reserve(page_limit);
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

    void sort_by_variant() {
        // TD: DO WE NEED SORT?
//      std::sort(data.begin(), data.end(), [](const VariantsPair& p1, const VariantsPair& p2) {
//        // Because of how the LDServer creates VariantPairs, it is guaranteed that the first variant in the pair
//        // is always before the second variant on the chromosome.
//        if (p1.variant1 == p2.variant1) {
//            return p1.position2 < p2.position2;
//        }
//        else {
//            return p1.position1 < p2.position1;
//        }
//      });
    }

    /**
     * Restrict this result object down to only variants given in the set.
     * @param variants Set of variants.
     */
    template<template <typename...> class C>
    void filter_by_variants(const C<string>& variants) {
        unsigned int idx = 0u, new_idx = 0u;
        unsigned int n_variants = this->variants.size();
        unordered_map<uint32_t, uint32_t> new_index_order;
        auto new_index_order_it = new_index_order.end();
        for (auto variants_it = this->variants.begin(); variants_it != this->variants.end(); ++idx) {
            if (variants.find(variants_it->name) == variants.end()) {
                variants_it = this->variants.erase(variants_it);
            } else {
                new_index_order.emplace(idx, new_idx++);
                ++variants_it;
            }
        }
        if (n_variants == this->variants.size()) { // nothing was remove i.e. no change
            return;
        }
        unordered_map<uint32_t, vector<Correlation>> new_correlations;
        unsigned int n_new_correlations = 0u;
        for (auto correlations_it = this->correlations.begin(); correlations_it != this->correlations.end(); ++correlations_it) {
            if ((new_index_order_it = new_index_order.find(correlations_it->first)) == new_index_order.end()) {
                continue;
            }
            new_idx = new_index_order_it->second;
            for (auto buddies_it = correlations_it->second.begin(); buddies_it != correlations_it->second.end(); ) {
                if ((new_index_order_it = new_index_order.find(buddies_it->variant_idx)) == new_index_order.end()) {
                    buddies_it = correlations_it->second.erase(buddies_it);
                } else {
                    buddies_it->variant_idx = new_index_order_it->second;
                    ++n_new_correlations;
                    ++buddies_it;
                }
            }
            if (correlations_it->second.empty()) {
                continue;
            }
            new_correlations.emplace(new_idx, move(correlations_it->second));
        }
        this->correlations = new_correlations;
        this->n_correlations = n_new_correlations;
//        vector<VariantsPair> new_data;
//        copy_if(data.begin(), data.end(), back_inserter(new_data), [&](const VariantsPair& p) -> bool {
//          return (variants.find(p.variant1) != variants.end()) && (variants.find(p.variant2) != variants.end());
//        });
//        data = new_data;
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
     * This function is called by LDQueryResult::get_json() to provide as part of the "next URL".
     *
     * @return
     */
    string get_last() const {
        if (has_next()) {
            return string(to_string(last_cell) + ":" + to_string(last_i) + ":" + to_string(last_j) + ":" + to_string(page));
        }
        return string("");
    }
    void clear_data() {
        variants.clear();
        index.clear();
        correlations.clear();
        n_correlations = 0u;
//        data.clear();
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
    string get_json(const string& url) {
        rapidjson::StringBuffer strbuf;
        rapidjson::Writer<rapidjson::StringBuffer> writer(strbuf);

        auto start = std::chrono::system_clock::now();
        writer.StartObject();
        writer.Key("data");
        writer.StartObject();
        writer.Key("variants");
        writer.StartArray();
        for (auto&& value: this->variants) {
            writer.String(value.name.c_str());
        }
        writer.EndArray();
        writer.Key("chromosomes");
        writer.StartArray();
        for (auto&& v: this->variants) {
            writer.String(v.chromosome.c_str());
        }
        writer.EndArray();
        writer.Key("positions");
        writer.StartArray();
        for (auto&& v: this->variants) {
            writer.Uint64(v.position);
        }
        writer.EndArray();
        writer.Key("correlation");

        writer.StartObject();
        for (auto&& c: correlations) {
            writer.Key(to_string(c.first).c_str());
            writer.StartArray();
            for (auto&& b: c.second) {
                writer.StartArray();
                writer.Uint(b.variant_idx);
                if (std::isnan(b.value)) {
                    writer.Null();
                } else {
                    writer.Double(b.value);
                }
                writer.EndArray();
            }
            writer.EndArray();
        }
        writer.EndObject();

//        writer.StartArray();
//        for (auto&& values : correlations) {
//            writer.StartArray();
//            for (auto&& value: values) {
//                if (std::isnan(value)) {
//                    writer.Null();
//                } else {
//                    writer.Double(value);
//                }
//            }
//            writer.EndArray();
//        }
//        writer.EndArray();
        writer.EndObject();
        writer.Key("error");
        writer.Null();
        writer.Key("next");
        if (is_last()) {
            writer.Null();
        } else {
            string link = url + "&last=" + get_last();
            writer.String(link.c_str());
        }
        writer.EndObject();
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double>  elapsed = end - start;
        std::cout << "Writing JSON to string elapsed time: " << elapsed.count() << " s\n";
        return strbuf.GetString();
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
