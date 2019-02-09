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
#include <boost/algorithm/string/predicate.hpp>

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

void Phenotypes::load_tab(const string &path, const ColumnTypeMap &types, size_t nrows) {
  ifstream input_file(path);
  string line;
  auto separator = regex("[ \t]");

  if (!input_file.good()) {
    throw std::invalid_argument("Cannot access file: " + path);
  }

  // Get header.
  vector<string> header;
  getline(input_file, line);
  copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(header));

  // Create vectors as necessary.
  ColumnType ct;
  string col;
  for (auto& kv : types) {
    col = kv.first;
    ct = kv.second;
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
  }

  // Convenience to lookup column type by index.
  vector<ColumnType> vec_types(types.size());
  for (int i = 0; i < header.size(); i++) {
    try {
      vec_types[i] = types.at(header[i]);
    }
    catch (...) {
      throw "Column " + header[i] + " did not exist in column type definitions";
    }
  }

  // Now read entire file, storing columns as we go.
  // We also now know the number of columns.
  vector<string> tokens;
  int i = 0;
  while (getline(input_file, line)) {
    // Skip commented lines.
    if (line.substr(0, 1) == "#") { continue; }
    if (line.empty()) { continue; }

    // Parse line.
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
    for (int j = 0; j < tokens.size(); j++) {
      col = header[j];
      switch (vec_types[j]) {
        case ColumnType::INTEGER:
          // For now store integers as doubles (as if they were floating point).
          // If we need space for some reason later, can make a separate storage for smaller ints.
        case ColumnType::FLOAT: {
          if (tokens[j] == "NA") {
            (*columns_float[col])(i) = arma::datum::nan;
          }
          else {
            (*columns_float[col])(i) = stod(tokens[j]);
          }
          break;
        }
        case ColumnType::TEXT: {
          columns_text[col]->emplace_back(tokens[j]);
          break;
        }
        case ColumnType::CATEGORICAL: {
          // If it's a NA value, we can skip parsing.
          if (tokens[j] == "NA") {
            (*columns_float[col])(i) = arma::datum::nan;
            continue;
          }

          // Have we seen this category before?
          auto it = map_level[col].find(tokens[j]);

          if (it != map_level[col].end()) {
            // We've seen this category before. Get its value.
            (*columns_float[col])(i) = it->second;
          }
          else {
            // New category level found.
            // Find the next category level.
            double cat_level = 0;

            // If the value is an integer, we'll just assume we should use their encoding.
            if (is_integer(tokens[j])) {
              cat_level = stod(tokens[j]);
            }
            else {
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
            map_level[col][tokens[j]] = cat_level;
            map_cat[col][cat_level] = tokens[j];
          }
          break;
        }
      }
    }

    tokens.clear();
    i++;
  }

  // In a tab-delimited file, we assume the first column contains the sample IDs.
  sample_ids = columns_text[header[0]];

  // Store a copy of the column types.
  column_types = types;
}

void Phenotypes::load_ped(const string &ped_path) {
  throw "Not yet implemented";
}

SharedVector<string> Phenotypes::as_text(const string &colname) {
  SharedVector<string> vec;
  switch (column_types[colname]) {
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

SharedArmaVec Phenotypes::as_vec(const string &colname) {
  SharedArmaVec vec;
  switch (column_types[colname]) {
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
    auto it = find(sample_ids->begin(), sample_ids->end(), sample);
    if (it != sample_ids->end()) {
      auto index = distance(sample_ids->begin(), it);
      indices.emplace_back(index);
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
  result->sigma2 = sigma2;
  result->pvalue = pvalue;

  return result;
}