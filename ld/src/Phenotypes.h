#ifndef LDSERVER_PHENOTYPES_H
#define LDSERVER_PHENOTYPES_H

#include <string>
#include <vector>
#include <armadillo>
#include "boost/variant.hpp"

using namespace std;

enum ColumnType : uint8_t {
  TEXT,
  CATEGORICAL,
  INTEGER,
  FLOAT
};

struct ScoreResult {
  double score_stat;
  double sigma2;
  double pvalue;
};

//template <typename T> T most_common(vector<T>& vec);
typedef boost::variant<vector<int64_t>, vector<double>, vector<string>> ColumnVariant;
typedef map<string, ColumnType> ColumnTypeMap;

template<class T> using SharedVector = shared_ptr<vector<T>>;
template<typename T> shared_ptr<vector<T>> make_shared_vector(vector<T>& v);

class Phenotypes {
private:
  map<string, shared_ptr<ColumnVariant>> columns;
  ColumnTypeMap column_types;
  SharedVector<string> sample_ids;
public:
  /**
   * Load a PED file and its associated DAT file.
   * Column types do not need to be provided, because the DAT file specifies them.
   * @param ped_path
   * @param dat_path
   */
  void load_ped(const string &ped_path, const string &dat_path);

  /**
   * Load a tab-delimited file.
   * Column types must be provided.
   * @param path
   * @param types
   */
  void load_tab(const string &path, const ColumnTypeMap &types);

  ColumnType get_column_type(const string& colname);
  shared_ptr<ColumnVariant> get_column(const string &colname);
  shared_ptr<arma::vec> as_vec(const string &colname);
  SharedVector<string> get_phenotypes();

  /**
   * Reduce this matrix of phenotypes to the provided list of samples (in the order given.)
   * Missing samples will have a NaN value for all phenotypes.
   * @param samples
   */
  void reorder(const vector<string> &samples);

  /**
   * Calculate score statistic, phenotypic variance, and p-value given a vector of genotypes.
   * @param genotypes
   * @param phenotype
   * @return
   */
  shared_ptr<ScoreResult> compute_score(arma::vec genotypes, const string &phenotype);
};

#endif //LDSERVER_PHENOTYPES_H