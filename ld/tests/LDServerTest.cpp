//#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <fstream>
#include <regex>
#include <map>
#include <boost/process/system.hpp>
#include "../src/LDServer.h"

using namespace std;

class LDServerTest: public::testing::Test {
protected:
    static boost::process::child redis_server;
    static redisContext* redis_cache;

    virtual ~LDServerTest() {}

    static void SetUpTestCase() {
        LDServerTest::redis_server = boost::process::child("../bin/redis-server --port 6379 --bind 127.0.0.1 --daemonize no --save \"\"");
        const char *hostname = "127.0.0.1";
        int port = 6379;
        this_thread::sleep_for(chrono::seconds(3));
        redis_cache = redisConnect(hostname, port);
    }

    static void TearDownTestCase() {
        if (redis_cache != nullptr) {
            redisFree(redis_cache);
        }
        LDServerTest::redis_server.terminate();
    }

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

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

    void load_samples(const string &file, vector<string>& samples) {
        ifstream input_file(file);
        string line;
        while (getline(input_file, line)) {
            samples.emplace_back(line);
        }
    }
};

boost::process::child LDServerTest::redis_server = boost::process::child();
redisContext* LDServerTest::redis_cache = nullptr;

TEST_F(LDServerTest, SAV_one_page) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(1000);

    server.set_file("chr22.test.sav");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
    }
}

TEST_F(LDServerTest, BCF_one_page) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(1000);

    server.set_file("chr22.test.bcf");
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
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
    ASSERT_TRUE(server.compute_region_ld("22", 51241101, 51241385, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
    }
}

TEST_F(LDServerTest, SAV_chrX_one_page) {
    map<string, double> goldstandard;
    vector<string> samples;
    this->load_region_goldstandard("region_ld_X_60100_60150.hap.ld", goldstandard);
    this->load_samples("EUR.samples.txt", samples);

    LDServer server;
    LDQueryResult result(1000);

    server.set_file("chrX.test.sav");
    ASSERT_TRUE(server.compute_region_ld("X", 60100, 60150, result));

    ASSERT_EQ(result.limit, 1000);
    ASSERT_EQ(result.get_last(), "");
    ASSERT_EQ(result.data.size(), goldstandard.size());
    for (auto&& entry : result.data) {
        string key(to_string(entry.position1) + "_" + to_string(entry.position2));
        ASSERT_NE(entry.variant1, "");
        ASSERT_NE(entry.variant2, "");
        ASSERT_EQ(goldstandard.count(key), 1);
        ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.000000001);
    }
}


TEST_F(LDServerTest, region_with_paging) {
    map<string, double> goldstandard;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(4);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_region_ld("22", 51241101, 51241385, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 4);
        ASSERT_LE(result.data.size(), 4);
        for (auto &&entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}


TEST_F(LDServerTest, variant_with_paging) {
    map<string, double> goldstandard;
    this->load_variant_goldstandard("variant_ld_22_51241101_vs_51241101_51241385.hap.ld", goldstandard);

    LDServer server;
    LDQueryResult result(2);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    while (server.compute_variant_ld("22:51241101_A/T", "22", 51241101, 51241385, result, LDServer::ALL_SAMPLES_KEY)) {
        ASSERT_LE(result.limit, 2);
        ASSERT_LE(result.data.size(), 2);
        for (auto&& entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.00000000001);
        }
        result_total_size += result.data.size();
    }
    ASSERT_EQ(result_total_size, goldstandard.size());
}

TEST_F(LDServerTest, AFR_region_with_paging) {
    map<string, double> goldstandard;
    vector<string> samples;
    this->load_region_goldstandard("region_ld_22_51241101_51241385.AFR.hap.ld", goldstandard);
    this->load_samples("AFR.samples.txt", samples);

    LDServer server;
    LDQueryResult result(2);
    int result_total_size = 0;

    server.set_file("chr22.test.sav");
    server.set_samples("AFR", samples);
    while (server.compute_region_ld("22", 51241101, 51241385, result, "AFR")) {
        ASSERT_LE(result.limit, 4);
        ASSERT_LE(result.data.size(), 4);
        for (auto &&entry : result.data) {
            string key(to_string(entry.position1) + "_" + to_string(entry.position2));
            ASSERT_NE(entry.variant1, "");
            ASSERT_NE(entry.variant2, "");
            ASSERT_EQ(goldstandard.count(key), 1);
            ASSERT_NEAR(goldstandard.find(key)->second , entry.rsquare, 0.0000000001);
        }
        result_total_size += result.data.size();
    }
}

TEST_F(LDServerTest, hiredis) {
    redisContext *context = nullptr;
    const char *hostname = "127.0.0.1";
    int port = 6379;
    struct timeval timeout = { 1, 500000 }; // 1.5 seconds
    context = redisConnectWithTimeout(hostname, port, timeout);
    ASSERT_NE(context, nullptr);
    ASSERT_EQ(context->err, 0);
    redisFree(context);
}

TEST_F(LDServerTest, segment_key) {
    uint64_t start_bp = 0u, stop_bp = 0u;

    Segment segment1("", 0, 1);
    ASSERT_EQ(16, segment1.get_key_size());
    memcpy(&start_bp, segment1.get_key(), 8);
    memcpy(&stop_bp, segment1.get_key() + 8, 8);
    ASSERT_EQ(start_bp, 0);
    ASSERT_EQ(stop_bp, 1);

    Segment segment2("chr22", 10, 20);
    ASSERT_EQ(21, segment2.get_key_size());
    string chromosome(segment2.get_key(), 5);
    ASSERT_EQ(chromosome, "chr22");
    memcpy(&start_bp, segment2.get_key() + 5 , 8);
    memcpy(&stop_bp, segment2.get_key() + 13, 8);
    ASSERT_EQ(start_bp, 10);
    ASSERT_EQ(stop_bp, 20);
}

TEST_F(LDServerTest, segment_cache) {
    RawSAV raw("chr22.test.sav");
    Segment segment1("22", 51241101, 51241385);
    raw.load(raw.get_samples(), segment1);
    segment1.save(redis_cache);

    Segment segment2("22", 51241101, 51241385);
    ASSERT_EQ(segment2.names.size(), 0);
    ASSERT_EQ(segment2.positions.size(), 0);
    segment2.load(redis_cache);

    ASSERT_THAT(segment1.names, ::testing::ContainerEq(segment2.names));
    ASSERT_THAT(segment1.positions, ::testing::ContainerEq(segment2.positions));
}

TEST_F(LDServerTest, cell_key) {
    uint64_t morton_code = 0;
    Cell cell1("", 3);
    ASSERT_EQ(8, cell1.get_key_size());
    memcpy(&morton_code, cell1.get_key(), 8);
    ASSERT_EQ(morton_code, 3);

    Cell cell2("chr22", 300);
    ASSERT_EQ(13, cell2.get_key_size());
    string chromosome(cell2.get_key(), 5);
    ASSERT_EQ(chromosome, "chr22");
    memcpy(&morton_code, cell2.get_key() + 5 , 8);
    ASSERT_EQ(morton_code, 300);
}

TEST_F(LDServerTest, cell_cache) {
    LDQueryResult result1(1000);
    LDQueryResult result2(1000);

    RawSAV raw("chr22.test.sav");
    Cell cell1("22", to_morton_code(512411, 512411));
    cell1.segment_i = make_shared<Segment>("22", 51241100, 51241199);
    raw.load(raw.get_samples(), *(cell1.segment_i));
    cell1.segment_i->save(redis_cache);
    cell1.compute();
    cell1.save(redis_cache);
    cell1.extract(51241100, 51241199, result1);

    Cell cell2("22", to_morton_code(512411, 512411));
    cell2.segment_i = make_shared<Segment>("22", 51241100, 51241199);
    cell2.segment_i->load(redis_cache);
    cell2.load(redis_cache);
    cell2.extract(51241100, 51241199, result2);

    ASSERT_THAT(cell1.segment_i->names, ::testing::ContainerEq(cell2.segment_i->names));
    ASSERT_THAT(cell1.segment_i->positions, ::testing::ContainerEq(cell2.segment_i->positions));
    ASSERT_EQ(result1.data.size(), result2.data.size());
    for (unsigned int i = 0; i < result1.data.size(); ++i) {
        ASSERT_EQ(result1.data[i].variant1, result2.data[i].variant1);
        ASSERT_EQ(result1.data[i].variant2, result2.data[i].variant2);
        ASSERT_EQ(result1.data[i].r, result2.data[i].r);
    }
}