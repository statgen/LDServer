#ifndef LDSERVER_PHENOTYPES_H
#define LDSERVER_PHENOTYPES_H

#include <utility>
#include <string>
#include <vector>
#include <armadillo>
#include "boost/variant.hpp"
#include "Types.h"
#include <stdexcept>

using namespace std;

//template <typename T> T most_common(vector<T>& vec);

class ColumnTypeMap {
protected:
  vector<pair<string, ColumnType>> types;
  map<string, ColumnType> ctmap;
public:
  inline void add(string name, ColumnType type) {
    types.emplace_back(make_pair(name, type));
    ctmap[name] = type;
  }
  inline ColumnType get_type(string name) { return ctmap.at(name); }
  inline auto size() const { return types.size(); }
  inline auto begin() const { return types.begin(); }
  inline auto end() const { return types.end(); }
};

class PhenotypeParseException : public std::runtime_error { using std::runtime_error::runtime_error; };

class Phenotypes {
protected:
  map<string, SharedArmaVec> columns_float;
  map<string, SharedVector<string>> columns_text;
  map<string, map<double, string>> map_cat;
  map<string, map<string, double>> map_level;

  string file_path;
  ColumnTypeMap column_types;
  SharedVector<string> sample_ids;
  map<string, uint64_t> sample_id_index;
public:
   /**
   * Load a phenotype file.
   * The file may be either:
   *   1. A tab-delimited file, with one phenotype per column. Must have a header row. File extension can be either
   *      .tab or .tab.gz.
   *   2. A PED-formatted file. Must have file extension .ped or .ped.gz. There must be an accompanying .dat or .dat.gz
   *      file in addition. You need only specify the path to the ped file.
   * @param path Path to tab or PED file.
   * @param types
   * @param nrows
   * @param delim File delimiter
   * @param sample_column Name of the column that contains the sample identifiers
   * @param analysis_cols Optional list of columns to actually parse (i.e. these columns will be used for analysis.)
   */
  void load_file(const string &path, const ColumnTypeMap &types, size_t nrows, const string& delim, const string& sample_column, SharedVector<string> analysis_cols = nullptr);

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
   * Get list of samples for which the phenotype has complete data.
   */
  shared_ptr<vector<string>> get_complete_samples(const string& phenotype);

  /**
   * Get number of non-missing observations for phenotype.
   * @param phenotype
   * @return
   */
  uint64_t get_nsamples(const string &phenotype);

  /**
   * Pretty print a summary of the current state of this object,
   * primarly for debugging purposes.
   */
  void pprint() const;
};

#endif //LDSERVER_PHENOTYPES_H