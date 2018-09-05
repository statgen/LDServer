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

    virtual ~LDServerTest() {}

    static void SetUpTestCase() {
        LDServerTest::redis_server = boost::process::child("../bin/redis-server --port 6379 --bind 127.0.0.1 --daemonize no --save \"\"");
    }

    static void TearDownTestCase() {
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
    redisContext *c = nullptr;
    const char *hostname = "127.0.0.1";
    int port = 6379;
    struct timeval timeout = { 1, 500000 }; // 1.5 seconds
    c = redisConnectWithTimeout(hostname, port, timeout);
    ASSERT_NE(c, nullptr);
    ASSERT_EQ(c->err, 0);
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

TEST_F(LDServerTest, segment_to_redis) {
    redisContext* context = nullptr;
    const char* hostname = "127.0.0.1";
    int port = 6379;
    struct timeval timeout = { 1, 500000 }; // 1.5 seconds
    context = redisConnectWithTimeout(hostname, port, timeout);
    ASSERT_NE(context, nullptr);
    ASSERT_EQ(context->err, 0);
    RawSAV raw("chr22.test.sav");
    Segment segment("22", 51241101, 51241385);
    raw.load(raw.get_samples(), segment);
    strstreambuf buffer;
    basic_ostream<char> os(&buffer);
    {
        cereal::BinaryOutputArchive oarchive(os);
        oarchive(segment);
    } // needs to exit the scope to flush the output - see cereal docs.
    redisReply* reply = nullptr;
    reply = (redisReply*)redisCommand(context, "SET %b %b", segment.get_key(), segment.get_key_size(), buffer.str(), buffer.pcount());
    ASSERT_NE(reply, nullptr);
    ASSERT_EQ(reply->type, REDIS_REPLY_STATUS);
    ASSERT_STREQ(reply->str, "OK");
    freeReplyObject(reply);
    redisFree(context);
}

TEST_F(LDServerTest, segment_from_redis) {
    redisContext* context = nullptr;
    const char* hostname = "127.0.0.1";
    int port = 6379;
    struct timeval timeout = { 1, 500000 }; // 1.5 seconds
    context = redisConnectWithTimeout(hostname, port, timeout);
    ASSERT_NE(context, nullptr);
    ASSERT_EQ(context->err, 0);

    RawSAV raw("chr22.test.sav");
    Segment segment_from_file("22", 51241101, 51241385);
    raw.load(raw.get_samples(), segment_from_file);

    Segment segment_from_redis("22", 51241101, 51241385);
    ASSERT_EQ(segment_from_redis.names.size(), 0);
    ASSERT_EQ(segment_from_redis.positions.size(), 0);

    redisReply* reply = nullptr;
    reply = (redisReply*)redisCommand(context, "GET %b", segment_from_redis.get_key(), segment_from_redis.get_key_size());
    ASSERT_EQ(reply->type, REDIS_REPLY_STRING);

    strstreambuf buffer(reply->str, reply->len);
    basic_istream<char> is(&buffer);
    {
        cereal::BinaryInputArchive iarchive(is);
        iarchive(segment_from_redis);
    }
    freeReplyObject(reply);

    ASSERT_THAT(segment_from_file.names, ::testing::ContainerEq(segment_from_redis.names));
    ASSERT_THAT(segment_from_file.positions, ::testing::ContainerEq(segment_from_redis.positions));
    redisFree(context);
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