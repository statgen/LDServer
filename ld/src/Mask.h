#ifndef LDSERVER_MASK_H
#define LDSERVER_MASK_H

#include <memory>
#include <string>
#include <map>
#include <vector>
#include <regex>
#include <iostream>
#include <set>
#include <algorithm>
#include <iterator>
#include <tabixpp.hpp>
#include "Types.h"

// TODO: move to Types.h
template<typename T>
struct VariantSort {
  bool operator()(const T& lhs, const T& rhs) const {
    return lhs.position < rhs.position;
  }
};

typedef std::set<VariantMeta, VariantSort<VariantMeta>> SortedVariantSet;

struct MaskGroup {
  std::string name;
  std::string chrom;
  uint64_t start;
  uint64_t stop;
  SortedVariantSet variants;

  std::shared_ptr<std::set<std::string>> get_variants() const;
};

class Mask {
  using group_iterator = std::map<std::string, MaskGroup>::const_iterator;

public:
  /**
   * Default constructor that loads the entire mask.
   * @param filepath
   */
  Mask(const std::string &filepath);

  /**
   * Constructor to load only a subset of the mask, using only regions of variants that overlap the given start/stop.
   * @param filepath
   * @param chrom
   * @param start
   * @param stop
   */
  Mask(const std::string &filepath, const std::string &chrom, uint64_t start, uint64_t stop);

  /**
   * Print out each group and its variants, mainly for debugging purposes.
   */
  void print_groups() const;

  /**
   * Get a list of group names.
   */
  std::shared_ptr<std::vector<std::string>> get_group_names() const;

  /**
   * Retrieve a group.
   */
  std::shared_ptr<MaskGroup> get_group(const std::string& group) const;

  /**
   * Functions for retrieving the variants in a group.
   * @param group The group name.
   */
  std::shared_ptr<std::set<std::string>> get_variant_set(const std::string &group) const;

  /**
   * Iterator functions.
   * The iterator will have first and second members, with the first member
   * being the name of the group, and the second member being a MaskGroup object.
   */
  group_iterator begin() const;
  group_iterator end() const;

protected:
  std::string name;
  uint64_t id;
  std::string description;

  /**
   * Store map from group name --> vector of variants.
   */
  std::map<std::string, MaskGroup> groups;

  /**
   * Loader for mask files. Expects the file to be tabixed and bgzipped.
   *
   * The first 4 columns should be:
   *   group name (for example, a gene)
   *   chromosome
   *   position of furthest upstream variant in the group (start)
   *   position of furthest downstream variant in the group (stop)
   *
   * The remainder of the line should be a tab-delimited list of variants in EPACTS format (chrom:pos_ref/alt).
   *
   * @param filepath
   * @param chrom
   * @param start
   * @param stop
   */
  void load_file(const std::string &filepath, const std::string &chrom = "", uint64_t start = 0, uint64_t stop = 0);
};

#endif //LDSERVER_MASK_H
