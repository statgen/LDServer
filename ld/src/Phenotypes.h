#ifndef LDSERVER_PHENOTYPES_H
#define LDSERVER_PHENOTYPES_H

#include <string>
#include <vector>
#include <armadillo>
#include "boost/variant.hpp"
#include "Types.h"

using namespace std;

//template <typename T> T most_common(vector<T>& vec);
typedef map<string, ColumnType> ColumnTypeMap;

class Phenotypes {
protected:
  map<string, SharedArmaVec> columns_float;
  map<string, SharedVector<string>> columns_text;
  map<string, map<double, string>> map_cat;
  map<string, map<string, double>> map_level;

  ColumnTypeMap column_types;
  SharedVector<string> sample_ids;
public:
  /**
   * Load a PED file and its associated DAT file.
   * The DAT file is assumed to be named as (basename of ped file).dat(.gz if ped was gzipped)
   * Column types do not need to be provided, because the DAT file specifies them.
   * @param ped_path
   */
  void load_ped(const string &ped_path);

  /**
   * Load a tab-delimited file.
   * Column types must be provided.
   * @param path
   * @param types
   */
  void load_tab(const string &path, const ColumnTypeMap &types, size_t nrows);

  SharedArmaVec as_vec(const string &colname);
  SharedVector<string> as_text(const string &colname);
  SharedVector<string> get_phenotypes();

  /**
   * Reduce this matrix of phenotypes to the provided list of samples (in the order given.)
   * Missing samples will have a NaN value for all phenotypes.
   * @param samples
   */
  void reorder(const vector<string> &samples);

  /**
   * Calculate score statistic, and p-value given a vector of genotypes.
   * @param genotypes
   * @param phenotype
   * @return
   */
  shared_ptr<ScoreResult> compute_score(arma::vec &genotypes, const string &phenotype);

  /**
   * Calculate phenotypic variance.
   * @param phenotype
   * @return
   */
  double compute_sigma2(const string &phenotype);

  /**
   * Get number of non-missing observations for phenotype.
   * @param phenotype
   * @return
   */
  uint64_t get_nsamples(const string &phenotype);
};

#endif //LDSERVER_PHENOTYPES_H