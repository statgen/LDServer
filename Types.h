#ifndef LDSERVER_TYPES_H
#define LDSERVER_TYPES_H

#include <string>
#include <vector>
#include <algorithm>
#include <regex>

using namespace std;

struct VariantMeta {
    string variant;
    string chromosome;
    string ref;
    string alt;
    uint64_t position;
    VariantMeta(const string& variant, const string& chromosome, const string& ref, const string& alt, uint64_t position):
            variant(variant), chromosome(chromosome), ref(ref), alt(alt), position(position) {}
    bool operator==(VariantMeta const& result) const { // needed by boost.python
        return variant.compare(result.variant) == 0;
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

struct VariantsPairLD {
    string variant1;
    string variant2;
    string chromosome1;
    string chromosome2;
    uint64_t position1;
    uint64_t position2;
    double r;
    double rsquare;
    VariantsPairLD() : variant1(""), variant2(""), chromosome1(""), chromosome2(""), position1(0ul), position2(0ul), r(0.0), rsquare(0.0) {}
    VariantsPairLD(const string& variant1, const string& chromosome1, uint64_t position1, const string& variant2, const string& chromosome2, uint64_t position2, double r, double rsquare):
            variant1(variant1), variant2(variant2), chromosome1(chromosome1), chromosome2(chromosome2), position1(position1), position2(position2), r(r), rsquare(rsquare) {}
    bool operator==(VariantsPairLD const& result) const { // needed by boost.python
        return (variant1.compare(result.variant1) == 0 && variant2.compare(result.variant2) == 0);
    }
};

struct LDQueryResult {
    uint32_t limit;
    uint64_t last_cell;
    int last_i;
    int last_j;
    int page;
    vector<VariantsPairLD> data;
    LDQueryResult(uint32_t page_limit): limit(page_limit), last_cell(0), last_i(-1), last_j(-1), page(0) {}
    LDQueryResult(const string& last) : limit(0), last_cell(0), last_i(-1), last_j(-1) {
        vector<std::string> tokens;
        copy(sregex_token_iterator(last.begin(), last.end(), regex(":"), -1), sregex_token_iterator(), back_inserter(tokens));
        if (tokens.size() > 3) {
            last_cell = std::stoull(tokens[0]);
            last_i = std::stoi(tokens[1]);
            last_j = std::stoi(tokens[2]);
            page = std::stoi(tokens[3]);
        }
    }
    bool has_next() const {
        return ((last_i >= 0) || (last_j >= 0));
    }
    bool is_last() const {
        return ((page > 0) && (last_i < 0) && (last_j < 0));
    }
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
};

#endif
