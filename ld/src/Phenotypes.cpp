#include "Phenotypes.h"
#include <vector>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;

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
template <typename T> T most_common(vector<T>& vec) {
  unordered_map<T, int> counts;
  int max = 0;
  T most_common;
  for (auto& v : vec) {
    int c = counts[v]++;
    if (c > max) {
      max = c;
      most_common = v;
    }
  }
  return most_common;
}

void Phenotypes::load_tab(const string &path) {
  ifstream input_file(path);
  string line;
  auto separator = regex("[ \t]");

  // Get header.
  vector<string> header;
  getline(input_file, line);
  copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(header));

  // Do a quick iteration of the first X lines to perform some basic column type inference.
  map<int, vector<ColumnType>> guesses;
  uint64_t line_counter = 0;
  while (getline(input_file, line)) {
    // Skip commented lines.
    if (line.substr(0,1) == "#") { continue; }
    if (line.empty()) { continue; }

    // Parse line.
    vector<string> tokens;
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));

    // Does the number of tokens match the number of elements in the header?
    if (header.size() != tokens.size()) {
      throw "Number of elements on row " + to_string(line_counter) + "did not match number in header";
    }

    for (int i = 0; i < tokens.size(); i++) {
      double tmp_double;
      int64_t tmp_int;

      // Try converting to floating point type first.
      try {
        tmp_double = stod(tokens.at(i));
      }
      catch (...) {
        // Can't convert to a double, so we have to leave it as a string.
        guesses[i].emplace_back(ColumnType::TEXT);
        continue;
      }

      // We at least know it's a number now. Could it be an integer?
      try {
        tmp_int = stol(tokens.at(i));
        if (tmp_double != tmp_int) {
          // Floating representation didn't match integer, so we can't safely store as integer.
          guesses[i].emplace_back(ColumnType::FLOAT);
        }
        else {
          guesses[i].emplace_back(ColumnType::INTEGER);
        }
      }
      catch (...) {
        guesses[i].emplace_back(ColumnType::FLOAT);
      }
    }

    if (line_counter > Phenotypes::N_LOOKAHEAD_LINES) {
      break;
    }

    line_counter++;
  }

  // Determine type of each column from our guesses, and create
  // our column store.
  ColumnType ct;
  string col;
  shared_ptr<VectorVariant> ptr;
  vector<ColumnType> vec_types(guesses.size());
  for (auto& kv : guesses) {
    ct = most_common(kv.second);
    col = header[kv.first];
    column_types[col] = ct;
    vec_types[kv.first] = ct;
    switch (ct) {
      case ColumnType::INTEGER:
        ptr = make_shared<VectorVariant>(vector<int64_t>());
        columns[col] = ptr;
        break;
      case ColumnType::FLOAT:
        ptr = make_shared<VectorVariant>(vector<double>());
        columns[col] = ptr;
        break;
      case ColumnType::TEXT:
        ptr = make_shared<VectorVariant>(vector<string>());
        columns[col] = ptr;
        break;
    }
  }

  // Now read entire file, storing columns as we go.
  // We also now know the number of columns.
  while (getline(input_file, line)) {
    // Skip commented lines.
    if (line.substr(0, 1) == "#") { continue; }
    if (line.empty()) { continue; }

    // Parse line.
    string col;
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
      }
    }
  }

  // In a tab-delimited file, we assume the first column contains the sample IDs.
  sample_ids = columns[header[0]];
}

void Phenotypes::load_ped(const string &ped_path, const string &dat_path) {
  throw "Not yet implemented";
}

shared_ptr<VectorVariant> Phenotypes::get_column(const string &colname) {
  return columns[colname];
}

ColumnType Phenotypes::get_column_type(const string& colname) {
  return column_types[colname];
}

shared_ptr<arma::vec> Phenotypes::as_float(const string &colname) {
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
  }

  return rvec;
}