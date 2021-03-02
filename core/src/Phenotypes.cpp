#include "Phenotypes.h"
#include <vector>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <armadillo>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/format.hpp>

using namespace std;

template<typename T> shared_ptr<vector<T>> make_shared_vector(vector<T>& v) { return make_shared<vector<T>>(v); }

inline bool is_integer(const std::string &s) {
  if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;

  char *p;
  strtol(s.c_str(), &p, 10);

  return (*p == 0);
}

/**
 * Find the most common element in a vector.
 * @tparam T
 * @param vec
 * @return
 */
//template <typename T> T most_common(vector<T>& vec) {
//  unordered_map<T, int> counts;
//  int max = 0;
//  T most_common;
//  for (auto& v : vec) {
//    int c = counts[v]++;
//    if (c > max) {
//      max = c;
//      most_common = v;
//    }
//  }
//  return most_common;
//}

void split(std::string const& original, const char& separator, std::vector<std::string>& results) {
  std::string::const_iterator start = original.begin();
  std::string::const_iterator end = original.end();
  std::string::const_iterator next = std::find(start, end, separator);
  while (next != end) {
    results.emplace_back(std::string(start, next));
    start = next + 1;
    next = std::find(start, end, separator);
  }
  results.emplace_back(std::string(start, next));
}

void Phenotypes::load_file(const string &path, const ColumnTypeMap &types, size_t nrows, const string& delim, const string& sample_column, SharedVector<string> analysis_cols) {
  ifstream input_file(path);
  string line;
  auto separator = regex(delim);

  if (!input_file.is_open()) {
    auto msg_full = boost::str(boost::format("Error opening phenotype file %s, error was: %s") % path % strerror(errno));
    auto msg_safe = string("Error opening phenotype file");
    throw LDServerGenericException(msg_safe).set_secret(msg_full);
  }

  this->file_path = path;

  // Do we need to read a header?
  bool is_ped = path.find(".ped") != string::npos;
  bool is_tab = !is_ped;
  if (is_tab) {
    // We need to read through the header line.
    getline(input_file, line);
  }

  // TODO: this can be cleaned up by just using the ColumnTypeMap directly and implementing the [] operator

  // If analysis_cols is empty, it means load all columns.
  if (analysis_cols == nullptr) {
    analysis_cols = make_shared<vector<string>>();
    transform(types.begin(), types.end(), back_inserter(*analysis_cols), [](const auto& p) { return p.first; });
  }
  else if (analysis_cols->empty()) {
    transform(types.begin(), types.end(), back_inserter(*analysis_cols), [](const auto& p) { return p.first; });
  }

  // Sample column must always be loaded
  if (find(analysis_cols->begin(), analysis_cols->end(), sample_column) == analysis_cols->end()) {
    analysis_cols->emplace_back(sample_column);
  }

    // Create vectors as necessary.
  vector<string> header;
  vector<uint64_t> parse_cols;
  ColumnTypeMap new_types;
  ColumnType ct;
  string col;
  uint64_t col_counter = -1;
  for (auto& kv : types) {
    col = kv.first;
    ct = kv.second;
    col_counter++;

    // We still want to know the name of every column in order in the file, even if we may not use them.
    header.emplace_back(col);

    if (find(analysis_cols->begin(), analysis_cols->end(), col) == analysis_cols->end()) {
      // Column is not meant to be analyzed
      continue;
    }

    switch (ct) {
      case ColumnType::INTEGER: {
        // Lazily store integer variables as floating point for now, can move to other storage later
        columns_float[col] = make_shared<arma::vec>(nrows);
        break;
      }
      case ColumnType::FLOAT: {
        columns_float[col] = make_shared<arma::vec>(nrows);
        break;
      }
      case ColumnType::TEXT: {
        columns_text[col] = make_shared<vector<string>>();
        break;
      }
      case ColumnType::CATEGORICAL: {
        // Same with integers, lazily store categorical variables encoded in floating point for now, move later
        columns_float[col] = make_shared<arma::vec>(nrows);
        break;
      }
    }

    new_types.add(col, ct);
    parse_cols.emplace_back(col_counter);
  }

  // Convenience to lookup column type by index.
  // Note this maps to every single column possible, not just the "for analysis" columns.
  vector<ColumnType> vec_types;
  transform(types.begin(), types.end(), back_inserter(vec_types), [](const auto& p) { return p.second; });

  // Now read entire file, storing columns as we go.
  // We also now know the number of columns.
  vector<string> tokens;
  int i = 0;
  while (getline(input_file, line)) {
    // Skip commented lines.
    if (line.substr(0, 1) == "#") { continue; }
    if (line.empty()) { continue; }

    // Parse line.
    split(line, delim[0], tokens);
    for (auto&& j : parse_cols) {
      col = header[j];
      string& val = tokens[j];

      try {
        switch (vec_types[j]) {
          case ColumnType::INTEGER:
            // For now store integers as doubles (as if they were floating point).
            // If we need space for some reason later, can make a separate storage for smaller ints.
          case ColumnType::FLOAT: {
            if ((val == "NA") || (val == ".") || (val == "")) {
              (*columns_float[col])(i) = arma::datum::nan;
            } else {
              (*columns_float[col])(i) = stod(val);
            }
            break;
          }
          case ColumnType::TEXT: {
            columns_text[col]->emplace_back(val);
            break;
          }
          case ColumnType::CATEGORICAL: {
            // If it's a NA value, we can skip parsing.
            if ((val == "NA") || (val == ".") || (val == "")) {
              (*columns_float[col])(i) = arma::datum::nan;
              continue;
            } else if (is_ped && ((val == "0") || (val == "-9"))) {
              // Binary traits in PED format specify 0 or -9 as missing value.
              (*columns_float[col])(i) = arma::datum::nan;
              continue;
            }

            // Have we seen this category before?
            auto it = map_level[col].find(val);

            if (it != map_level[col].end()) {
              // We've seen this category before. Get its value.
              (*columns_float[col])(i) = it->second;
            } else {
              // New category level found.
              // Find the next category level.
              double cat_level = 0;

              if (is_integer(val)) {
                if (is_ped) {
                  auto tmp = stoull(val);
                  if (tmp > 2) {
                    throw std::range_error(
                      "Categorical variables in PED files are expected to be coded 0=missing,1=unaffected,2=affected, found value: " +
                      val);
                  } else if (tmp == 1) {
                    cat_level = 0;
                  } else if (tmp == 2) {
                    cat_level = 1;
                  }
                } else {
                  // If the value is an integer, we'll just assume we should use their encoding.
                  cat_level = stod(val);
                }
              } else {
                if (!map_level[col].empty()) {
                  cat_level = (*max_element(
                    map_level[col].begin(),
                    map_level[col].end(),
                    [](const auto &v1, const auto &v2) {
                      return v1.second < v2.second;
                    }
                  )).second + 1;
                }
              }

              // Assign next level.
              (*columns_float[col])(i) = cat_level;

              // Remember level encoding
              map_level[col][val] = cat_level;
              map_cat[col][cat_level] = val;
            }
            break;
          }
        }
      }
      catch (std::exception& e) {
        auto fmt_msg = boost::format("Error reading phenotype file on line %i, column %i (%s), invalid value: %s");
        auto s_msg = boost::str(fmt_msg % i % j % header[j] % val);

        auto fmt_secret = boost::format("Bad phenotype file: %s, original exception was %s: %s");
        auto s_secret = boost::str(fmt_secret % path % typeid(e).name() % e.what());

        throw LDServerGenericException(s_msg).set_secret(s_secret);
      }
    }

    tokens.clear();
    i++;
  }

  if (input_file.bad()) {
    auto msg_full = boost::str(boost::format("Error reading phenotype file %s, error was: %s") % path % strerror(errno));
    auto msg_safe = string("Error while reading phenotype file");
    throw LDServerGenericException(msg_safe).set_secret(msg_full);
  }

  // In a tab-delimited file, we assume the first column contains the sample IDs.
  sample_ids = columns_text[sample_column];
  for (uint64_t i = 0; i < sample_ids->size(); i++) {
    sample_id_index[sample_ids->at(i)] = i;
  }

  // Store a copy of the column types.
  column_types = new_types;
}

SharedVector<string> Phenotypes::as_text(const string &colname) {
  SharedVector<string> vec;
  auto type = column_types.get_type(colname);
  switch (type) {
    case ColumnType::INTEGER:
    case ColumnType::FLOAT:
    case ColumnType::CATEGORICAL:
      throw "Cannot convert column " + colname + "to text";
    case ColumnType::TEXT:
      vec = columns_text[colname];
      break;
  }

  return vec;
}

// shared_ptr<arma::vec>
SharedArmaVec Phenotypes::as_vec(const string &colname) {
  SharedArmaVec vec;
  auto type = column_types.get_type(colname);
  switch (type) {
    case ColumnType::INTEGER:
    case ColumnType::FLOAT:
    case ColumnType::CATEGORICAL:
      vec = columns_float[colname];
      break;
    case ColumnType::TEXT:
      throw "Cannot convert text column to floating point";
  }

  return vec;
}

/**
 * Utility function to help reordering an arbitrary container given an ordering of indices into it.
 * @tparam T
 * @param container
 * @param indices Ordered list of indices into container. An index of -1 means insert NaN.
 * @return shared_ptr to newly created / reordered container
 */
template<typename T> shared_ptr<T> reorder_container(const T& container, const vector<int64_t>& indices) {
  auto new_container = make_shared<T>(indices.size());

  for (int i = 0; i < indices.size(); i++) {
    const int64_t &ind = indices[i];
    if (ind == -1) {
      (*new_container)[i] = arma::datum::nan;
    }
    else {
      (*new_container)[i] = container[ind];
    }
  }

  return new_container;
}

void Phenotypes::reorder(const vector<string> &samples) {
  // Find indices of existing samples
  vector<int64_t> indices;
  for (auto& sample : samples) {
    // Find index of requested sample in our original list of samples.
    auto it = sample_id_index.find(sample);
    if (it != sample_id_index.end()) {
      indices.emplace_back(it->second);
    }
    else {
      // Store a sentinel value to denote this sample didn't exist.
      indices.emplace_back(-1);
    }
  }

  // Reorder samples accordingly. For samples that didn't exist, we need to store a NaN in the proper location.
  for (auto& kv : column_types) {
    switch (kv.second) {
      case ColumnType::INTEGER:
      case ColumnType::FLOAT: {
        columns_float[kv.first] = reorder_container(*columns_float[kv.first], indices);
        break;
      }
      case ColumnType::TEXT: {
        columns_text[kv.first] = reorder_container(*columns_text[kv.first], indices);
        break;
      }
      case ColumnType::CATEGORICAL: {
        columns_float[kv.first] = reorder_container(*columns_float[kv.first], indices);
        break;
      }
    }
  }

  // Set samples
  vector<string> copy_samples = samples;
  this->sample_ids = make_shared_vector<string>(copy_samples);
  for (uint64_t i = 0; i < sample_ids->size(); i++) {
    sample_id_index[sample_ids->at(i)] = i;
  }
}

SharedVector<string> Phenotypes::get_phenotypes() {
  vector<string> phenos;
  for (auto& kv : column_types) {
    phenos.emplace_back(kv.first);
  }
  return make_shared_vector<string>(phenos);
}

/**
 * Calculate score statistic and other relevant statistics.
 * Assumes that the matrix has already been reordered according to the genotypes' samples.
 * @param genotypes
 * @param phenotype
 * @return shared_ptr to ScoreResult
 */
shared_ptr<ScoreResult> Phenotypes::compute_score(arma::vec &genotypes, const string &phenotype) {
  // First get vector of phenotype.
  auto &pheno_vec = *as_vec(phenotype);

  // Find the mean of the genotypes in order to center them
  double mean = arma::mean(genotypes);

  // Create mean-centered genotype vector
  arma::vec geno_vec(genotypes - mean);

  // Calculate score statistic.
  double score_stat = arma::dot(geno_vec, pheno_vec);

  // Calculate sigma2.
  // The second parameter (1) tells arma to divide by N, not N-1.
  double sigma2 = arma::var(pheno_vec, 1);

  // Calculate p-value.
  double t_stat = score_stat / sqrt(sigma2);
  double pvalue = 2 * arma::normcdf(-fabs(t_stat));

  // Create return object.
  auto result = make_shared<ScoreResult>();
  result->score_stat = score_stat;
  result->pvalue = pvalue;

  return result;
}

/**
 * Get list of samples for which the phenotype has complete data.
 */
shared_ptr<vector<string>> Phenotypes::get_complete_samples(const string& phenotype) {
  // Figure out which phenotype rows have non-missing data.
  auto pheno_vec = as_vec(phenotype);
  arma::uvec index_nonmiss = find_finite(*pheno_vec);
  auto samples = make_shared<vector<string>>();

  for (uint64_t i = 0; i < index_nonmiss.n_elem; i++) {
    samples->emplace_back((*sample_ids)[index_nonmiss[i]]);
  }

  return samples;
}

double Phenotypes::compute_sigma2(const string &phenotype) {
  auto &pheno_vec = *as_vec(phenotype);
  arma::vec nonmiss_pheno = pheno_vec.elem(find_finite(pheno_vec));
  return arma::var(nonmiss_pheno, 1);
}

uint64_t Phenotypes::get_nsamples(const string &phenotype) {
  auto &pheno_vec = *as_vec(phenotype);
  arma::uvec non_finite = arma::find_nonfinite(pheno_vec);
  return pheno_vec.n_elem - non_finite.n_elem;
}

void Phenotypes::pprint() const {
  cout << "Loaded file: " << this->file_path << endl;
  cout << "Number of columns: " << this->column_types.size() << endl;
  cout << "Column types: " << endl;
  uint64_t lim = 5;
  for (auto&& p : this->column_types) {
    auto colname = p.first;
    auto ctype = p.second;
    cout << "- " + colname + ": " + to_string(ctype) << endl;

    switch (ctype) {
      case ColumnType::INTEGER:
      case ColumnType::CATEGORICAL:
      case ColumnType::FLOAT: {
        auto vec = *columns_float.at(colname);
        arma::uvec non_finite = arma::find_nonfinite(vec);
        cout << "  " << "Number of elements: " << to_string(vec.n_elem) << endl;
        cout << "  " << "Number of non-missing elements: " << to_string(vec.n_elem - non_finite.n_elem) << endl;

        cout << "  " << "First few elements:";
        uint64_t n_elem = vec.n_elem;
        uint64_t enum_n = std::min(n_elem, lim);
        for (int i = 0; i < enum_n; i++) {
          cout << " " << to_string(vec[i]);
        }
        cout << endl;

        break;
      }
      case ColumnType::TEXT: {
        auto vec = *columns_text.at(colname);
        cout << "  " << "Number of elements: " << to_string(vec.size()) << endl;

        cout << "  " << "First few elements:";
        uint64_t n_elem = vec.size();
        uint64_t enum_n = std::min(n_elem, lim);
        for (int i = 0; i < enum_n; i++) {
          cout << " " << vec[i];
        }
        cout << endl;

        break;
      }
    }

    cout << endl;
  }
}