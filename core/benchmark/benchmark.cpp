#include <benchmark/benchmark.h>
#include "../src/LDServer.h"
#include "../src/ScoreCovarianceRunner.h"
#include "../src/SummaryStatisticsLoader.h"
#include "../src/RaremetalSummaryStatisticsLoader.h"
#include "../src/MetastaarSummaryStatisticsLoader.h"

static void BM_LDSERVER_REGIONLD_RSQUARE(benchmark::State& state) {
  LDServer server(100);
  LDQueryResult result(100000);
  server.set_file("chr22.test.sav");

  for (auto _ : state) {
    // Compute
    server.compute_region_ld("22", 50244251, 51244237, correlation::LD_RSQUARE, result);
  }
}

BENCHMARK(BM_LDSERVER_REGIONLD_RSQUARE)->Unit(benchmark::kMillisecond)->Iterations(3);

static void BM_LDSERVER_JSON_CLASSIC(benchmark::State& state) {
  LDServer server(100);
  LDQueryResult result(100000);
  server.set_file("chr22.test.sav");

  // Compute
  server.compute_region_ld("22", 50244251, 51244237, correlation::LD_RSQUARE, result);

  for (auto _ : state) {
    // Get JSON
    auto json = result.get_json_classic("blah");

    // Parse back out JSON
    rapidjson::Document doc;
    doc.Parse(json.c_str());
  }
}

BENCHMARK(BM_LDSERVER_JSON_CLASSIC)->Unit(benchmark::kMillisecond)->Iterations(3);

static void BM_LDSERVER_JSON_COMPACT(benchmark::State& state) {
  LDServer server(100);
  LDQueryResult result(100000);
  server.set_file("chr22.test.sav");

  // Compute
  server.compute_region_ld("22", 50244251, 51244237, correlation::LD_RSQUARE, result);

  for (auto _ : state) {
    // Get JSON
    auto json = result.get_json_compact("blah");

    // Parse back out JSON
    rapidjson::Document doc;
    doc.Parse(json.c_str());
  }
}

BENCHMARK(BM_LDSERVER_JSON_COMPACT)->Unit(benchmark::kMillisecond)->Iterations(3);

static void BM_METASTAAR_LOADER(benchmark::State& state) {
  for (auto _ : state) {
    MetastaarSummaryStatisticsLoader loader(
      {
        "test.qt.segment1.metastaar.sumstat.parquet",
        "test.qt.segment2.metastaar.sumstat.parquet"
      },
      {
        "test.qt.segment1.metastaar.cov.parquet",
        "test.qt.segment2.metastaar.cov.parquet"
      }
    );

    loader.load_region("1", 4957, 5143);
  }
}

BENCHMARK(BM_METASTAAR_LOADER)->Unit(benchmark::kMillisecond)->Iterations(3);

static void BM_RAREMETAL_LOADER(benchmark::State& state) {
  for (auto _ : state) {
    RaremetalSummaryStatisticsLoader loader({"test.smallchunk.MetaScore.assoc.gz"}, {"test.smallchunk.MetaCov.assoc.gz"});
    loader.load_region("1", 2, 307);
  }
}

BENCHMARK(BM_RAREMETAL_LOADER)->Unit(benchmark::kMillisecond)->Iterations(3);

BENCHMARK_MAIN();