#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <regex>
#include <map>
#include "../LDServer.h"

using namespace std;

class LDServerTest: public::testing::Test {
protected:
    virtual ~LDServerTest() {}
    virtual void SetUp() {}
    virtual void TearDown() {}

    void load_region_goldstandard(const string &file, map<string, double> &values) {
        ifstream input_file(file);
        string line;
        vector<string> tokens;
        getline(input_file, line); // skip header;
        while (getline(input_file, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), regex("[ \t]"), -1), sregex_token_iterator(), back_inserter(tokens));
            values.emplace(tokens.at(1) + "_" + tokens.at(2), stod(tokens.at(4)));
            tokens.clear();
        }
    }

    void load_variant_goldstandard(const string &file, map<string, double> &values) {
        ifstream input_file(file);
        string line;
        vector<string> tokens;
        getline(input_file, line); // skip header;
        while (getline(input_file, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), regex("[ \t]"), -1), sregex_token_iterator(), back_inserter(tokens));
            values.emplace(tokens.at(1) + "_" + tokens.at(3), stod(tokens.at(5)));
            tokens.clear();
        }
    }
};


TEST_F(LDServerTest, SAV_one_page) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(1000);

    server.set_file("chr22.test.sav");
    server.compute_region_ld("22", 51241101, 51241385, result);

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
    }
}


TEST_F(LDServerTest, VCF_one_page) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(1000);

    server.set_file("chr22.test.vcf.gz");
    server.compute_region_ld("22", 51241101, 51241385, result);

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
    }
}


TEST_F(LDServerTest, region_with_paging) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(4);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    do {
        server.compute_region_ld("22", 51241101, 51241385, result, LDServer::ALL_SAMPLES_KEY);
        ASSERT_LE(result.limit, 4);
        ASSERT_LE(result.data.size(), 4);
        for (auto &&entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
        }
        result_total_size += result.data.size();
    } while (result.has_next());
    ASSERT_EQ(result_total_size, goldstandard.size());
}


TEST_F(LDServerTest, variant_with_paging) {
    map<string, double> goldstandard;
    this->load_variant_goldstandard("variant_ld_22_51241101_vs_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(2);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    do {
        server.compute_variant_ld("22:51241101_A/T", "22", 51241101, 51241385, result, LDServer::ALL_SAMPLES_KEY);
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
        for (auto&& entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
        }
        result_total_size += result.data.size();
    } while (result.has_next());
    ASSERT_EQ(result_total_size, goldstandard.size());
}