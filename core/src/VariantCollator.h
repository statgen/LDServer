#ifndef LDSERVER_VARIANTCOLLATOR_H
#define LDSERVER_VARIANTCOLLATOR_H

#include <vector>
#include <unordered_map>
#include "Types.h"
#include "Raw.h"
#include <tabixpp.hpp>
#include "RaremetalSummaryStatisticsLoader.h"
#include "MetastaarSummaryStatisticsLoader.h"
#include <savvy/reader.hpp>
#include <savvy/sav_reader.hpp>
#include <savvy/vcf_reader.hpp>
#include <savvy/region.hpp>
#include <savvy/site_info.hpp>

enum class VariantFileFormat {VCF, SAVVY, RAREMETAL, METASTAAR};

std::string to_string(VariantFileFormat& fmt);

class VariantCollator {
protected:
  VariantFileFormat format;
  std::unordered_map<std::string, std::string> chrom_file;
  std::map<std::string, MetastaarFileIntervalTree> score_tree;
  std::map<std::string, MetastaarFileIntervalTree> cov_tree;

  std::shared_ptr<std::vector<VariantMeta>> variants;

  void read_variants_genotype_file(std::string chrom, uint64_t start, uint64_t end);
  void read_variants_raremetal_file(std::string chrom, uint64_t start, uint64_t end);
  void read_variants_metastaar_file(std::string chrom, uint64_t start, uint64_t end);

public:
  VariantCollator(std::vector<std::string> genotype_files, VariantFileFormat format);
  VariantCollator(std::vector<std::string> score_files, std::vector<std::string> cov_files, VariantFileFormat format);
  std::shared_ptr<std::vector<VariantMeta>> get_variants(std::string chrom, uint64_t start, uint64_t end);
};


#endif //LDSERVER_VARIANTCOLLATOR_H
