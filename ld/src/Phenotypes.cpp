#include "Phenotypes.h"
#include <vector>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;

template<typename T> shared_ptr<vector<T>> make_shared_vector(vector<T>& v) { return make_shared<vector<T>>(v); }

//inline bool isInteger(const std::string & s) {
//  if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;
//
//  char *p;
//  strtol(s.c_str(), &p, 10);
//
//  return (*p == 0);
//}

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

void Phenotypes::load_tab(const string &path, const ColumnTypeMap &types) {
  ifstream input_file(path);
  string line;
  auto separator = regex("[ \t]");

  // Get header.
  vector<string> header;
  getline(input_file, line);
  copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(header));

  // Allocate storage from column types.
  ColumnType ct;
  string col;
  shared_ptr<ColumnVariant> ptr;
  for (auto& kv : types) {
    col = kv.first;
    ct = kv.second;
    switch (ct) {
      case ColumnType::INTEGER:
        ptr = make_shared<ColumnVariant>(vector<int64_t>());
        columns[col] = ptr;
        break;
      case ColumnType::FLOAT:
        ptr = make_shared<ColumnVariant>(vector<double>());
        columns[col] = ptr;
        break;
      case ColumnType::TEXT:
        ptr = make_shared<ColumnVariant>(vector<string>());
        columns[col] = ptr;
        break;
      case ColumnType::CATEGORICAL:
        ptr = make_shared<ColumnVariant>(vector<string>());
        columns[col] = ptr;
        break;
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
  while (getline(input_file, line)) {
    // Skip commented lines.
    if (line.substr(0, 1) == "#") { continue; }
    if (line.empty()) { continue; }

    // Parse line.
    vector<string> tokens;
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
    for (int i = 0; i < tokens.size(); i++) {
      col = header[i];
      switch (vec_types[i]) {
        case ColumnType::INTEGER:
          boost::get<vector<int64_t>>(*columns[col]).emplace_back(stol(tokens[i]));
          break;
        case ColumnType::FLOAT:
          boost::get<vector<double>>(*columns[col]).emplace_back(stod(tokens[i]));
          break;
        case ColumnType::TEXT:
          boost::get<vector<string>>(*columns[col]).emplace_back(tokens[i]);
          break;
        case ColumnType::CATEGORICAL:
          boost::get<vector<string>>(*columns[col]).emplace_back(tokens[i]);
          break;
      }
    }
  }

  // In a tab-delimited file, we assume the first column contains the sample IDs.
  sample_ids = make_shared<vector<string>>(boost::get<vector<string>>(*columns[header[0]]));

  // Store a copy of the column types.
  column_types = types;
}

void Phenotypes::load_ped(const string &ped_path, const string &dat_path) {
  throw "Not yet implemented";
}

shared_ptr<ColumnVariant> Phenotypes::get_column(const string &colname) {
  return columns[colname];
}

ColumnType Phenotypes::get_column_type(const string& colname) {
  return column_types[colname];
}

shared_ptr<arma::vec> Phenotypes::as_vec(const string &colname) {
  shared_ptr<arma::vec> rvec;
  arma::vec v;
  switch (get_column_type(colname)) {
    case ColumnType::INTEGER:
      {
        auto &ref = boost::get<vector<int64_t>>(*columns[colname]);
        vector<double> tmp(ref.begin(), ref.end());
        v = arma::vec(tmp);
        rvec = make_shared<arma::vec>(v);
      }
      break;
    case ColumnType::FLOAT:
      {
        v = arma::vec(boost::get<vector<double>>(*columns[colname]));
        rvec = make_shared<arma::vec>(v);
      }
      break;
    case ColumnType::TEXT:
      throw "Cannot convert text column to floating point";
    case ColumnType::CATEGORICAL:
      throw "Cannot convert categorical column to floating point";
  }

  return rvec;
}

void Phenotypes::reorder(const vector<string> &samples) {
  size_t n_new = samples.size();

  // Find indices of existing samples
  vector<uint64_t> indices;
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
  for (auto& kv : columns) {
    decltype(kv.second) new_col;

    for (auto i : indices) {

    }

  }

  // Set samples
  vector<string> copy_samples = samples;
  this->sample_ids = make_shared_vector<string>(copy_samples);
}

SharedVector<string> Phenotypes::get_phenotypes() {
  vector<string> phenos;
  for (auto& kv : columns) {
    phenos.emplace_back(kv.first);
  }
  return make_shared_vector<string>(phenos);
}

shared_ptr<ScoreResult> Phenotypes::compute_score(arma::vec genotypes, const string& phenotype) {

}