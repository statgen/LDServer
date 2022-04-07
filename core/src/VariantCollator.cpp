#include "VariantCollator.h"
using namespace std;

VariantCollator::VariantCollator(vector<string> genotype_files, VariantFileFormat format) {
  this->format = format;

  if (!(format == VariantFileFormat::SAVVY || format == VariantFileFormat::VCF)) {
    throw LDServerGenericException(boost::str(boost::format("Unknown file format requested when collating list of variants: %s") % to_string(format)));
  }

  // Create our own copy of the genotype file paths.
  chrom_file.clear();
  score_tree.clear();

  // Figure out which file contains which chromosome.
  for (auto& file : genotype_files) {
    shared_ptr<Raw> raw = RawFactory::create(file);
    for (auto& chrom: raw->get_chromosomes()) {
      this->chrom_file.emplace(chrom, file);
    }
  }
}

VariantCollator::VariantCollator(std::vector<std::string> score_files, std::vector<std::string> cov_files, VariantFileFormat format) {
  // Create our own copy of the genotype file paths.
  chrom_file.clear();
  score_tree.clear();
  cov_tree.clear();

  this->format = format;
  if (format == VariantFileFormat::RAREMETAL) {
    for (auto& f : score_files) {
      Tabix tb(const_cast<string&>(f));
      for (auto& chrom : tb.chroms) {
        chrom_file[chrom] = f;
      }
    }
  }
  else if (format == VariantFileFormat::METASTAAR) {
    typedef Interval<uint64_t, MetastaarParquetMetadata> MetaInterval;

    // Create map from chromosome to list of genomic intervals covered by each score statistic file.
    // Each MetaInterval contains the start/end position of the interval covered by the file, in addition to a copy of the
    // parquet metadata for each file.
    map<string, vector<MetaInterval>> score_map;
    for (auto& f : score_files) {
      auto meta = read_parquet_metadata(f);
      MetaInterval intv(meta.region_start, meta.region_mid, meta);
      score_map[meta.chrom].emplace_back(intv);
    }

    // Construct an interval tree for each chromosome representing the score statistic file intervals.
    for (const auto& p : score_map) {
      score_tree.emplace(make_pair(p.first, p.second));
    }

    // Same as above, only now created for each covariance matrix file.
    map<string, vector<MetaInterval>> cov_map;
    for (auto& f : cov_files) {
      auto meta = read_parquet_metadata(f);
      MetaInterval intv(meta.region_start, meta.region_mid, meta);
      cov_map[meta.chrom].emplace_back(intv);
    }

    for (const auto& p : cov_map) {
      cov_tree.emplace(make_pair(p.first, p.second));
    }
  }
  else {
    throw LDServerGenericException(boost::str(boost::format("Unknown file format requested when collating list of variants: %s") % to_string(format)));
  }
}

void VariantCollator::read_variants_metastaar_file(std::string chrom, uint64_t start, uint64_t end) {
  variants = make_shared<vector<VariantMeta>>();

  // Pull out interval tree for particular chromosome.
  const MetastaarFileIntervalTree* chrom_score_tree;
  const MetastaarFileIntervalTree* chrom_cov_tree;

  try {
    chrom_score_tree = &score_tree.at(chrom);
    chrom_cov_tree = &cov_tree.at(chrom);
  }
  catch (std::out_of_range& e) {
    throw LDServerGenericException(
      boost::str(boost::format("Chromosome %s not present in score stat files") % chrom)
    );
  }

  auto sorter = [](auto &i1, auto& i2) { return i1.start < i2.start; };

  // Retrieve list of score stat files needed as well.
  auto score_overlaps = chrom_score_tree->findOverlapping(start, end);
  if (score_overlaps.empty()) {
    throw LDServerGenericException(boost::str(boost::format("Region %s:%s-%s did not overlap any MetaSTAAR summary stat (score) file") % chrom % start % end));
  }
  std::sort(score_overlaps.begin(), score_overlaps.end(), sorter);

  for (uint64_t block = 0; block < score_overlaps.size(); block++) {
    auto score_int = score_overlaps[block];

    // For each file found, extract data.
    // Note: we need only load data where MAF > 0 and MAF > cov_maf_cutoff. Only variants matching those two criteria
    // are stored in the MetaSTAAR covariance files.
    std::shared_ptr<arrow::io::ReadableFile> score_infile;
    PARQUET_ASSIGN_OR_THROW(score_infile, arrow::io::ReadableFile::Open(score_int.value.filepath));
    parquet::StreamReader score_reader{parquet::ParquetFileReader::Open(score_infile)};
    int ncols = score_reader.num_columns();
    //int64_t nrows = score_reader.num_rows();
    int ncovariates = ncols - 10;

    #ifndef NDEBUG
    cout << boost::str(boost::format("Loading score statistics from file %s") % score_int.value.filepath) << endl;
    #endif

    // Get the covariance file corresponding to this same range.
    const auto cov_overlaps = chrom_cov_tree->findOverlapping(score_int.start, score_int.stop);
    if (cov_overlaps.size() > 1) {
      throw LDServerGenericException("Multiple MetaSTAAR covariance files overlapped a region covered by one score statistic file, should be one-to-one mapping")
        .set_secret(boost::str(boost::format("Score stat file was '%s' and region %s:%s-%s") % score_int.value.filepath % chrom % score_int.start % score_int.stop));
    }

    // MetaSTAAR cov files only store variants with MAF > 0 and MAF < maf_cutoff.
    // The summary stat / score file is a superset of those variants.
    const double& maf_cutoff = cov_overlaps[0].value.cov_maf_cutoff;

    // Temporary storage while reading each row of parquet file. We unfortunately absolutely need to pull each
    // value while reading or it will cause a parquet reader exception.
    string row_chrom;
    uint32_t pos;
    string ref;
    string alt;
    uint32_t alt_AC;
    uint32_t MAC;
    double maf;
    uint32_t N;
    double U;
    double V;

    double tmp;

    while (!score_reader.eof()) {
      // Note: you must read all values, not reading all values and sending EndRow will result in an exception
      score_reader >> row_chrom >> pos >> ref >> alt >> alt_AC >> MAC >> maf >> N >> U >> V;
      for (int i = 0; i < ncovariates; i++) {
        score_reader >> tmp;
      }
      score_reader >> parquet::EndRow;

      if ((maf < 0) || (maf >= maf_cutoff)) {
        // Variants with MAF < 0 or >= cutoff are not stored in the GtG file (the "cov" file.)
        continue;
      }

      if (pos >= start) {
        if (pos > end) {
          break;
        }

        variants->emplace_back(row_chrom, ref, alt, pos);
      }
    }
  }
}

void VariantCollator::read_variants_raremetal_file(std::string chrom, uint64_t start, uint64_t end) {
  variants = make_shared<vector<VariantMeta>>();

  string score_path;
  try {
    score_path = chrom_file.at(chrom);
  }
  catch (std::out_of_range& e) {
    throw NoVariantsInRange(
      boost::str(boost::format("Chromosome %s not present in score statistics files") % chrom)
    );
  }

  auto detected_format = detectScoreCovFormat(score_path);

#ifndef NDEBUG
  cout << boost::str(boost::format("Loading score statistics on chromosome %s from file %s") % chrom % score_path) << endl;
#endif

  Tabix tbfile(const_cast<string&>(score_path));
  string region = chrom + ":" + to_string(start) + "-" + to_string(end);

  bool has_chrom = find(tbfile.chroms.begin(), tbfile.chroms.end(), chrom) != tbfile.chroms.end();
  if (!has_chrom) {
    throw NoVariantsInRange("Chromosome " + chrom + " not found within score statistic file");
  }

  string line;
  vector<string> tokens;

  if ((!chrom.empty()) && (start != 0) && (end != 0)) {
    tbfile.setRegion(region);
  }

  const ScoreColumnSpec* cols;
  if (detected_format == ScoreCovFormat::RVTEST) {
    cols = &SCORE_COLUMNS_RVTEST;
  }
  else if (detected_format == ScoreCovFormat::RAREMETAL) {
    cols = &SCORE_COLUMNS_RAREMETAL;
  }

  uint64_t scores_read = 0;
  while (tbfile.getNextLine(line)) {
    auto separator = regex("[ \t]");
    copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));

    // Extract information from line
    try {
      string lineno = boost::str(boost::format("<lineno: %d>") % scores_read);
      variants->emplace_back(
        tokens.at(cols->colChrom),
        tokens.at(cols->colRef),
        tokens.at(cols->colAlt),
        extract_numeric<unsigned long>(spstoul, tokens.at(cols->colPos), cols->colPos, score_path, lineno)
      );
    }
    catch (LDServerGenericException& e) {
      throw;
    }
    catch(...) {
      throw LDServerGenericException(
        "Invalid value detected while parsing score statistic file"
      ).set_secret(
        boost::str(boost::format("File was: %s, offending line was:\n %s") % score_path % line)
      );
    }

    scores_read++;
    tokens.clear();
  }

  if (scores_read == 0) {
    throw NoVariantsInRange(
      boost::str(boost::format("No score statistics loaded within genomic region %s:%i-%i") % chrom % start % end)
    );
  }
}

void VariantCollator::read_variants_genotype_file(string chrom, uint64_t start, uint64_t end) {
  string filepath;
  try {
    filepath = chrom_file.at(chrom);
  }
  catch (std::out_of_range& e) {
    throw LDServerGenericException(
      boost::str(boost::format("Chromosome %s not present in genotype files") % chrom)
    );
  }

  bool is_vcf = filepath.find(".vcf") != string::npos;
  bool is_bcf = filepath.find(".bcf") != string::npos;
  bool is_sav = filepath.find(".sav") != string::npos;

  if (!(is_vcf || is_bcf || is_sav)) {
    throw std::invalid_argument("File " + filepath + " has unsupported format");
  }

  savvy::indexed_reader reader(filepath, {chrom, start, end}, savvy::fmt::gt);
  savvy::variant<vector<float>> var;

  // TODO: there is currently no way (due to compression) to seek through the genotypes and only read the variant site info
  // so the genotypes end up being loaded from disk, which wastes some time. It's possible we could separately index the
  // variants into a separate parquet file or something of that sort during server startup but may not be worth the effort.
  while (reader >> var) {
    variants->emplace_back(var.chromosome(), var.ref(), var.alt(), var.position());
  }
}

shared_ptr<vector<VariantMeta>> VariantCollator::get_variants(string chrom, uint64_t start, uint64_t end) {
  variants = make_shared<vector<VariantMeta>>();

  bool fmt_genotypes = (this->format == VariantFileFormat::VCF) || (this->format == VariantFileFormat::SAVVY);
  if (fmt_genotypes) {
    read_variants_genotype_file(chrom, start, end);
  }
  else if (this->format == VariantFileFormat::RAREMETAL) {
    read_variants_raremetal_file(chrom, start, end);
  }
  else if (this->format == VariantFileFormat::METASTAAR) {
    read_variants_metastaar_file(chrom, start, end);
  }

  return variants;
}