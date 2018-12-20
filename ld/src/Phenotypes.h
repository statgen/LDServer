#ifndef LDSERVER_PHENOTYPES_H
#define LDSERVER_PHENOTYPES_H

#include <string>
#include <vector>
#include <armadillo>
#include "boost/variant.hpp"

using namespace std;

enum class ColumnType {
  TEXT,
  FLOAT,
  INTEGER,
};

template <typename T> T most_common(vector<T>& vec);
typedef boost::variant<vector<int64_t>, vector<double>, vector<string>> VectorVariant;

class Phenotypes {
private:
  map<string, shared_ptr<VectorVariant>> columns;
  map<string, ColumnType> column_types;
  shared_ptr<VectorVariant> sample_ids;
  static const uint64_t N_LOOKAHEAD_LINES = 100;
public:
  // Loading methods
  void load_ped(const string &ped_path, const string &dat_path);
  void load_tab(const string &path);
  ColumnType get_column_type(const string& colname);
  shared_ptr<VectorVariant> get_column(const string &colname);
  shared_ptr<arma::vec> as_float(const string &colname);
};

#endif //LDSERVER_PHENOTYPES_H