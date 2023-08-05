#ifndef LDSERVER_METASTAARSUMMARYSTATISTICSLOADER_H
#define LDSERVER_METASTAARSUMMARYSTATISTICSLOADER_H

#include "SummaryStatisticsLoader.h"

/*
 Struct to represent the metadata stored within each MetaSTAAR parquet file.

 There are two parquet files per segment of the genome: a single variant statistic file, which contains information per
 variant such as chrom/pos/ref/alt, MAC/MAF, and score statistic. Appended at the end of the primary columns are extra
 columns representing the GtU matrix (where U is the half hat matrix; there is 1 column per covariate in the null model.)

 MetaSTAAR chunks the genome up into segments of some fixed size, called the segment size. Each segment has a start,
 mid, and end position (sometimes referred to as "region" start, mid, end.)

 Each segment though may or may not contain genetic variants, depending on the particular dataset. There may be no known
 variants there, or the genotyping array did not have variants typed there, etc.

 The score files look like:
                                                                                               GtU
                                                                                       ┌────────────────┐
                                                                                       ▼                ▼
                       chr  pos ref alt  alt_AC  MAC      MAF    N         U         V        1         2
    Segment 1 start ──►  1    1   C   A      53   53 0.010566 2508 -3.615483 12.264961 0.515152 -0.060252
                         1    2   A   G      57   57 0.011364 2508 -4.244432 13.173521 0.553539 -0.002675
    Segment 1 mid   ──►  1    3   A   C      52   52 0.010367 2508  1.089355 12.032375 0.504133  0.104684
                       ──────────────────────────────────────────────────────────────────────────────────
    Segment 2 start ──►  1    4   G   T      66   66 0.013158 2508 -3.910195 15.667785 0.641294 -0.047769
                         1    5   G   C      84   84 0.016746 2508  3.675083 21.090993 0.815441  0.033904
    Segment 2 mid   ──►  1    6   G   C      56   56 0.011164 2508  0.155301 13.420570 0.543676  0.016430

 Note that each segment's score file only contains the variants from the start to midpoint of each segment, and not the end.

 The covariance matrices look like:

       Segment   Segment
        Start    End
            │    │
            ▼    ▼
            123456
           ┌──────┐
          1│xxx...│
          2│.xxx..│
Segment──►3│..xxx.│
  Mid      └──────┘
            123456
              ▲
              │
           Segment
             Mid

 Note that each covariance matrix is rectangular, in order to store sliding windows of <segment size> covariances for
 each row variant. The rows only extend to the region midpoint, but the columns extend to the region end. In order to
 lookup information for variants past the midpoint, the next segment's score file needs to be loaded.
 */
struct MetastaarParquetMetadata {
  std::string filepath;
  std::string chrom;
  uint64_t region_start = 0;
  uint64_t region_mid = 0;
  uint64_t region_end = 0;
  uint64_t nrows = 0;
  uint64_t ncols = 0;
  double cov_maf_cutoff = 0;
};

MetastaarParquetMetadata read_parquet_metadata(const std::string& s);

using MetastaarFileIntervalTree = IntervalTree<uint64_t, MetastaarParquetMetadata>;

/**
 * Loader for MetaSTAAR summary statistic files.
 *
 * MetaSTAAR separates the final covariance matrix into
 */
class MetastaarSummaryStatisticsLoader : public SummaryStatisticsLoader {
protected:
  // Maps from chromosome -> interval tree of (start pos, end pos) for MetaSTAAR segmented files.
  std::map<std::string, MetastaarFileIntervalTree> score_tree;
  std::map<std::string, MetastaarFileIntervalTree> cov_tree;

  shared_ptr<LDQueryResult> cov_result;
  shared_ptr<ScoreStatQueryResult> score_result;
  uint64_t nsamples;
public:
  MetastaarSummaryStatisticsLoader(const std::vector<std::string>& score_vec, const std::vector<std::string>& cov_vec);
  void load_region(const std::string& chromosome, uint64_t start, uint64_t stop);
  shared_ptr<LDQueryResult> getCovResult();
  shared_ptr<ScoreStatQueryResult> getScoreResult();
  double getSigma2();
  uint64_t getNumSamples();
};

#endif
