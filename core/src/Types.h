#ifndef LDSERVER_TYPES_H
#define LDSERVER_TYPES_H

#include <string>
#include <vector>
#include <algorithm>
#include <regex>
#include <codecvt>
#include <locale>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <set>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/writer.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <boost/format.hpp>
#include <armadillo>

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

enum class FilterValueType {STRING, DOUBLE};

/**
 * Class to represent a variant filter.
 *
 * This would be much easier with templates, but that makes the C++/python bindings much more difficult (at least
 * with boost python. It's possible pybind11 improves the situation.)
 *
 * Instead, here we opt for an ersatz RTTI/union type system that works easily with python, but makes code complexity
 * worse for the C++ parts that must use it.
 */
struct VariantFilter {
  std::string op;
  std::string field;

  FilterValueType type;
  std::string value_string;
  double value_double = numeric_limits<double>::quiet_NaN();

  void set_value(const std::string& v) {
    value_string = v;
    type = FilterValueType::STRING;
  }

  void set_value(const double& v) {
    value_double = v;
    type = FilterValueType::DOUBLE;
  }

  template<typename T>
  T get_value() {
    switch (type) {
      case FilterValueType::STRING:
        return value_string;
      case FilterValueType::DOUBLE:
        return value_double;
    }
  }

  bool operator==(const VariantFilter& other) {
    bool both_nan = std::isnan(value_double) && std::isnan(other.value_double);
    bool one_nan = std::isnan(value_double) ^ std::isnan(other.value_double);

    if (one_nan) { return false; }
    if (!both_nan) {
      if (value_double != other.value_double) {
        return false;
      }
    }

    return (op == other.op) &&
           (field == other.field) &&
           (value_string == other.value_string);
  }
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
    vector<VariantsPair> data;
    LDQueryResult(uint32_t page_limit): limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0) {
        data.reserve(page_limit);
    }

    /**
     * Construct an LDQueryResult
     * @param page_limit
     * @param last This appears to be a string of the form "last_cell:last_i:last_j:page". last_cell is the morton code
     *   of the last cell that was retrieved.
     */
    LDQueryResult(uint32_t page_limit, const string& last) : limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0) {
        data.reserve(page_limit);
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
      std::sort(data.begin(), data.end(), [](const VariantsPair& p1, const VariantsPair& p2) {
        // Because of how the LDServer creates VariantPairs, it is guaranteed that the first variant in the pair
        // is always before the second variant on the chromosome.
        if (p1.variant1 == p2.variant1) {
            return p1.position2 < p2.position2;
        }
        else {
            return p1.position1 < p2.position1;
        }
      });
    }

    /**
     * Restrict this result object down to only variants given in the set.
     * @param variants Set of variants.
     */
    template<template <typename...> class C>
    void filter_by_variants(const C<string>& variants) {
        vector<VariantsPair> new_data;
        copy_if(data.begin(), data.end(), back_inserter(new_data), [&](const VariantsPair& p) -> bool {
          return (variants.find(p.variant1) != variants.end()) && (variants.find(p.variant2) != variants.end());
        });
        data = new_data;
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
        data.clear();
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
        rapidjson::Document document;
        document.SetObject();
        rapidjson::Document::AllocatorType& allocator = document.GetAllocator();

        rapidjson::Value data(rapidjson::kObjectType);

        rapidjson::Value variant1(rapidjson::kArrayType);
        rapidjson::Value chromosome1(rapidjson::kArrayType);
        rapidjson::Value position1(rapidjson::kArrayType);
        rapidjson::Value variant2(rapidjson::kArrayType);
        rapidjson::Value chromosome2(rapidjson::kArrayType);
        rapidjson::Value position2(rapidjson::kArrayType);
        rapidjson::Value correlation(rapidjson::kArrayType);

        for (auto&& p: this->data) {
            variant1.PushBack(rapidjson::StringRef(p.variant1.c_str()), allocator);
            chromosome1.PushBack(rapidjson::StringRef(p.chromosome1.c_str()), allocator);
            position1.PushBack(p.position1, allocator);
            variant2.PushBack(rapidjson::StringRef(p.variant2.c_str()), allocator);
            chromosome2.PushBack(rapidjson::StringRef(p.chromosome2.c_str()), allocator);
            position2.PushBack(p.position2, allocator);
            if (std::isnan(p.value)) { // nan is not allowed by JSON, so we replace it with null
                correlation.PushBack(rapidjson::Value(), allocator);
            } else {
                correlation.PushBack(p.value, allocator);
            }
        }

        data.AddMember("variant1", variant1, allocator);
        data.AddMember("chromosome1", chromosome1, allocator);
        data.AddMember("position1", position1, allocator);
        data.AddMember("variant2", variant2, allocator);
        data.AddMember("chromosome2", chromosome2, allocator);
        data.AddMember("position2", position2, allocator);
        data.AddMember("correlation", correlation, allocator);

        document.AddMember("data", data, allocator);
        document.AddMember("error", rapidjson::Value(), allocator);
        if (is_last()) {
            document.AddMember("next", rapidjson::Value(), allocator);
        } else {
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
  double sigma2 = numeric_limits<double>::quiet_NaN(); // assume all score stats calculated against the same phenotype
  double nsamples = numeric_limits<double>::quiet_NaN();

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

  std::shared_ptr<std::set<std::string>> get_variants() const {
    auto vs = make_shared<std::set<std::string>>();
    transform(
      data.begin(),
      data.end(),
      inserter(*vs, vs->begin()),
      [](const auto& v) {
        return v.variant;
      }
    );
    return vs;
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

  /**
   * Restrict this object down to only variants passing a specific filter.
   * @param f VariantFilter object specifying the field (maf, pvalue, etc.), operator (gte is >=, lte is <=), and the value
   *    to compare against.
   */
  void filter(const VariantFilter& f) {
    vector<ScoreResult> new_data;
    for (auto& it : data) {
      if (f.field == "maf") {
        double maf = min(it.alt_freq, 1-it.alt_freq);

        if (f.op == "gte") {
          if (maf >= f.value_double) {
            new_data.emplace_back(it);
          }
        }

        if (f.op == "lte") {
          if (maf <= f.value_double) {
            new_data.emplace_back(it);
          }
        }
      }

      else if (f.field == "pvalue") {
        if (f.op == "gte") {
          if (it.pvalue >= f.value_double) {
            new_data.emplace_back(it);
          }
        }

        if (f.op == "lte") {
          if (it.pvalue <= f.value_double) {
            new_data.emplace_back(it);
          }
        }
      }
    }

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
