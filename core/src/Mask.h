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
#include <boost/format.hpp>
#include "Types.h"

// TODO: move to Types.h
template<typename T>
struct VariantSort {
  bool operator()(const T& lhs, const T& rhs) const {
    return lhs.position < rhs.position;
  }
};

// TODO: could we just refer to these types by string, instead of having to define each one as a hardcoded enum?

// Enum for type of variant group (gene, region, other?)
enum VariantGroupType : uint8_t {
  GENE,
  REGION
};

// Enum for type of identifier for each group
enum GroupIdentifierType : uint8_t {
  ENSEMBL
};

typedef std::set<VariantMeta, VariantSort<VariantMeta>> SortedVariantSet;

struct VariantGroup {
  std::string name;
  std::string chrom;
  uint64_t start;
  uint64_t stop;
  SortedVariantSet variants;

  std::shared_ptr<std::set<std::string>> get_variants() const;
};

class Mask {
  using group_iterator = std::map<std::string, VariantGroup>::const_iterator;

public:
  /**
   * Default constructor that loads the entire mask.
   * @param filepath
   */
  Mask(const std::string &filepath, const uint64_t id, VariantGroupType group_type, GroupIdentifierType ident_type);

  /**
   * Constructor to load only a subset of the mask, using only regions of variants that overlap the given start/stop.
   * @param filepath
   * @param chrom
   * @param start
   * @param stop
   */
  Mask(const std::string &filepath, const uint64_t id, VariantGroupType group_type, GroupIdentifierType ident_type, const std::string &chrom, uint64_t start, uint64_t stop);

  /**
   * Print out each group and its variants, mainly for debugging purposes.
   * @param group_limit Limit the number of groups printed in total.
   * @param variant_limit Limit the number of variants printed out for each group.
   */
  void print_groups(const uint64_t& group_limit, const uint64_t& variant_limit) const;

  /**
   * Get a list of group names.
   */
  std::shared_ptr<std::vector<std::string>> get_group_names() const;

  /**
   * Retrieve a group.
   */
  std::shared_ptr<VariantGroup> get_group(const std::string& group) const;

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

  /**
   * Getters/setters
   */
  inline uint64_t get_id() const { return id; };

  inline VariantGroupType get_group_type() const { return group_type; };
  inline void set_group_type(VariantGroupType group_type) { this->group_type = group_type; }

  inline GroupIdentifierType get_identifier_type() const { return identifier_type; };
  inline void set_identifier_type(GroupIdentifierType identifier_type) { this->identifier_type = identifier_type; }

  /**
   * Operators
   */
  bool operator==(const Mask& other) const;

protected:
  uint64_t id;                         // Unique string for this mask
  std::string description;                // Text description of what variant filters this mask represents
  VariantGroupType group_type;            // The type of group (is it a gene, a region, etc?)
  GroupIdentifierType identifier_type;    // The identifier type of the group (ENSEMBL ID, REFSEQ?)

  /**
   * Store map from group name --> vector of variants.
   */
  std::map<std::string, VariantGroup> groups;

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
