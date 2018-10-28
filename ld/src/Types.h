#ifndef LDSERVER_TYPES_H
#define LDSERVER_TYPES_H

#include <string>
#include <vector>
#include <algorithm>
#include <regex>
#include <codecvt>
#include <locale>
#include <cmath>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/writer.h>
#include <cereal/external/rapidjson/stringbuffer.h>

using namespace std;

enum correlation : uint8_t {
    LD_R,
    LD_RSQUARE,
    COV
};

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

#endif
