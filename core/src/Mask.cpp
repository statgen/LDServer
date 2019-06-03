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

void VariantGroup::add_variant(const std::string& variant) {
  VariantMeta vm(variant);
  variants.emplace(vm);
  this->chrom = vm.chromosome;
  this->start = this->start == -1 ? vm.position : std::min(this->start, vm.position);
  this->stop = this->stop == -1 ? vm.position : std::max(this->stop, vm.position);
}

void Mask::load_file(const string &filepath, const string &chrom, uint64_t start, uint64_t stop) {
  if (start <= 0) { throw std::invalid_argument("Mask starting position was < 0"); }
  if (stop <= 0) { throw std::invalid_argument("Mask stop position was < 0"); }

  Tabix tbfile(const_cast<string&>(filepath));
  string region = chrom + ":" + to_string(start) + "-" + to_string(stop);

  bool has_chrom = find(tbfile.chroms.begin(), tbfile.chroms.end(), chrom) != tbfile.chroms.end();
  if (!has_chrom) {
    throw std::range_error("Chromosome " + chrom + " not found within mask file");
  }

  string line;
  vector<string> tokens;

  if ((!chrom.empty()) && (start != 0) && (stop != 0)) {
    tbfile.setRegion(region);
  }

  uint64_t groups_added = 0;
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
    groups_added++;

    tokens.clear();
  }

  if (groups_added == 0) {
    throw std::range_error(
      boost::str(boost::format("No groups loaded within genomic region %s:%i-%i for mask %s") % chrom % start % stop % id)
    );
  }
}

Mask::Mask(const string& filepath, const uint64_t id, VariantGroupType group_type, GroupIdentifierType ident_type) {
  this->id = id;
  this->group_type = group_type;
  this->identifier_type = ident_type;
  load_file(filepath);
}

Mask::Mask(const string &filepath, const uint64_t id, VariantGroupType group_type, GroupIdentifierType ident_type, const string &chrom, uint64_t start, uint64_t stop) {
  this->id = id;
  this->group_type = group_type;
  this->identifier_type = ident_type;
  load_file(filepath, chrom, start, stop);
}

void Mask::print_groups(const uint64_t& group_limit, const uint64_t& variant_limit) const {
  uint64_t g = 0;
  for (auto&& kv : groups) {
    cout << kv.second.name << endl;
    cout << "Chrom: " << kv.second.chrom << endl;
    cout << "Start: " << kv.second.start << endl;
    cout << "Stop: " << kv.second.stop << endl;
    cout << "Variants: " << endl;
    uint64_t v = 0;
    for (auto&& vmeta : kv.second.variants) {
      cout << "  " + vmeta.variant << endl;
      if (v > variant_limit) {
        break;
      }
      v++;
    }

    if (g > group_limit) {
      break;
    }
    g++;
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

bool Mask::operator==(const Mask& other) const {
  return id == other.id;
}