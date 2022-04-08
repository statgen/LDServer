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

/**
 * Enum for various variant file formats. Variant site information can be stored in VCF/BCF, SAVVY, RAREMETAL (also used
 * for rvtest), and METASTAAR files.
 */
enum class VariantFileFormat {VCF, SAVVY, RAREMETAL, METASTAAR};

/**
 * Convert an object of enum class VariantFileFormat into a string description.
 * @param fmt - VariantFileFormat object
 * @return - std::string description
 */
std::string to_string(VariantFileFormat& fmt);

/**
 * Class to retrieve lists of variants from various filetypes (vcf, savvy, raremetal/rvtest, metastaar.)
 */
class VariantCollator {
protected:
  VariantFileFormat format;

  // Used to keep track of chromosome -> score statistic file (RAREMETAL or rvtest only) or genotype file.
  std::unordered_map<std::string, std::string> chrom_file;

  // Data structures used to keep track of MetaSTAAR score and covariance files.
  std::map<std::string, MetastaarFileIntervalTree> score_tree;
  std::map<std::string, MetastaarFileIntervalTree> cov_tree;

  std::shared_ptr<std::vector<VariantMeta>> variants;

  void read_variants_genotype_file(std::string chrom, uint64_t start, uint64_t end);
  void read_variants_raremetal_file(std::string chrom, uint64_t start, uint64_t end);
  void read_variants_metastaar_file(std::string chrom, uint64_t start, uint64_t end);

public:
  /**
   * Constructor to use when given a list of genotype files (VCF/BCF, SAVVY), one file per chromosome.
   * Each file must have a tabix index on disk (chr1.vcf.gz & chr1.vcf.gz.tbi must exist). Only the genotype file though
   * needs to be provided in the genotype_files vector.
   * @param genotype_files - vector of genotype file paths
   * @param format - format of the file. Use VariantFileFormat enum class.
   */
  VariantCollator(std::vector<std::string> genotype_files, VariantFileFormat format);

  /**
   * Constructor to use when given a set of summary statistic files (score statistic files & covariance matrix files.)
   * The covariance files are unfortunately required for MetaSTAAR, as the files contain important
   * metadata needed to filter the list of variants in the score statistic file down to the proper subset (those with MAF cutoff
   * below what was used to generate the covariance matrix.)
   * @param score_files - vector of score statistic file paths
   * @param cov_files - vector of covariance file paths
   * @param format - format of the file (use VariantFileFormat enum class.) If given rvtest file, use VariantFileFormat::RAREMETAL.
   */
  VariantCollator(std::vector<std::string> score_files, std::vector<std::string> cov_files, VariantFileFormat format);

  /**
   * Get a list of variants in a region.
   * @param chrom
   * @param start
   * @param end
   * @return Vector of VariantMeta objects, one for each variant within the requested region.
   */
  std::shared_ptr<std::vector<VariantMeta>> get_variants(std::string chrom, uint64_t start, uint64_t end);
};


#endif //LDSERVER_VARIANTCOLLATOR_H
