#include "Mask.h"
using namespace std;

shared_ptr<set<string>> VariantGroup::get_variants() const {
  auto vs = make_shared<set<string>>();
  transform(
    variants.begin(),
    variants.end(),
    inserter(*vs, vs->begin()),
    [](const VariantMeta& vm) {
      return vm.variant;
    }
  );
  return vs;
}

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

    // Create group object
    VariantGroup group;
    group.name = tokens[0];
    group.chrom = tokens[1];
    group.start = stoull(tokens[2]);
    group.stop = stoull(tokens[3]);

    // Extract variants
    SortedVariantSet vset;
    transform(
      tokens.begin() + 4,
      tokens.end(),
      inserter(vset, vset.begin()),
      [](const string &str) {
        return VariantMeta(str);
      }
    );
    group.variants = vset;

    // Store group to map
    groups.emplace(make_pair(tokens[0], group));

    tokens.clear();
  }
}

Mask::Mask(const string& filepath, const std::string& name, VariantGroupType group_type) {
  this->name = name;
  this->group_type = group_type;
  load_file(filepath);
}

Mask::Mask(const string &filepath, const std::string& name, VariantGroupType group_type, const string &chrom, uint64_t start, uint64_t stop) {
  this->name = name;
  this->group_type = group_type;
  load_file(filepath, chrom, start, stop);
}

void Mask::print_groups() const {
  for (auto&& kv : groups) {
    cout << kv.second.name << endl;
    cout << "Chrom: " << kv.second.chrom << endl;
    cout << "Start: " << kv.second.start << endl;
    cout << "Stop: " << kv.second.stop << endl;
    cout << "Variants: " << endl;
    for (auto&& vmeta : kv.second.variants) {
      cout << "  " + vmeta.variant << endl;
    }
  }
}

shared_ptr<set<string>> Mask::get_variant_set(const string &group) const {
  auto iter = groups.find(group);
  if (iter == groups.end()) {
    throw out_of_range("Group " + group + "not found in mask file");
  }

  auto ptr = make_shared<set<string>>();
  transform(
    iter->second.variants.begin(),
    iter->second.variants.end(),
    inserter(*ptr, (*ptr).begin()),
    [](const VariantMeta& vm) {
      return vm.variant;
    }
  );

  return ptr;
}

shared_ptr<vector<string>> Mask::get_group_names() const {
  auto group_names = make_shared<vector<string>>(groups.size());
  for (auto&& kv : groups) {
    group_names->emplace_back(kv.first);
  }
  return group_names;
}

shared_ptr<VariantGroup> Mask::get_group(const string& group) const {
  auto iter = groups.find(group);
  if (iter == groups.end()) {
    throw out_of_range("Group " + group + "not found in mask file");
  }

  auto ptr = make_shared<VariantGroup>(iter->second);
  return ptr;
}

Mask::group_iterator Mask::begin() const {
  return groups.begin();
}

Mask::group_iterator Mask::end() const {
  return groups.end();
}