#ifndef LDSERVER_MASK_H
#define LDSERVER_MASK_H

#include <string>
#include <map>
#include <vector>
#include <regex>
#include <iostream>
#include <set>
#include <tabixpp.hpp>

class Mask {
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
  Mask(const std::string &filepath, const std::string &chrom, uint64_t &start, uint64_t &stop);

  /**
   * Print out each group and its variants, mainly for debugging purposes.
   */
  void print_groups() const;

  /**
   * Functions for retrieving the variants in a group.
   * @param group The group name.
   */
  shared_ptr<std::set<std::string>> get_variant_set(const std::string &group) const;
  shared_ptr<std::vector<std::string>> get_variant_vector(const std::string &group) const;

protected:
  /**
   * Store map from group name --> vector of variants.
   */
  std::map<std::string, std::vector<std::string>> groups;

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
