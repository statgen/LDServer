#include "Mask.h"
using namespace std;

void Mask::load_file(const string &filepath, const string &chrom, uint64_t start, uint64_t stop) {
  Tabix tbfile(const_cast<string&>(filepath));
  string region = chrom + ":" + to_string(start) + "-" + to_string(stop);

  string line;
  vector<string> tokens;

  if ((!chrom.empty()) && (start != 0) && (stop != 0)) {
    tbfile.setRegion(region);
  }


  while (tbfile.getNextLine(line)) {
    auto separator = regex("[ \t]");
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));

    // Store list of variants by group name.
    vector<string> variants(tokens.begin() + 4, tokens.end());
    groups.emplace(make_pair(tokens[0], variants));

    tokens.clear();
  }
}

Mask::Mask(const string& filepath) {
  load_file(filepath);
}

Mask::Mask(const string &filepath, const string &chrom, uint64_t &start, uint64_t &stop) {
  load_file(filepath, chrom, start, stop);
}

void Mask::print_groups() const {
  for (auto&& kv : groups) {
    cout << kv.first << endl;
    for (auto&& variant : kv.second) {
      cout << "  " + variant << endl;
    }
  }
}

shared_ptr<set<string>> Mask::get_variant_set(const std::string &group) const {
  auto variant_iter = groups.find(group);
  if (variant_iter == groups.end()) {
    throw out_of_range("Group " + group + "not found in mask file");
  }

  auto& variants = variant_iter->second;
  auto ptr = make_shared<set<string>>(variants.begin(), variants.end());
  return ptr;
}

shared_ptr<vector<string>> Mask::get_variant_vector(const std::string &group) const {
  auto variant_iter = groups.find(group);
  if (variant_iter == groups.end()) {
    throw out_of_range("Group " + group + "not found in mask file");
  }

  auto ptr = make_shared<vector<string>>(variant_iter->second);
  return ptr;
}