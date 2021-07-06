#include "ScoreCovarianceRunner.h"
#include <sys/stat.h>
using namespace std;
using namespace rapidjson;

const uint32_t MAX_UINT32 = numeric_limits<uint32_t>::max();
const uint32_t INITIAL_RESULT_SIZE = 10000000;

bool is_file(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

void ScoreCovarianceConfig::pprint() const {
  cout << boost::format("Region: %s:%i-%i") % chrom % start % stop << endl;
  cout << "Genotype dataset ID: " << genotype_dataset_id << endl;
  cout << "Genotype files: " << endl;
  if (genotype_files.empty()) {
    cout << ".. ** NO GENOTYPE FILES FOUND **" << endl;
  }
  else {
    for (auto &&f : genotype_files) {
      cout << ".. " << f << endl;
    }
  }
  cout << "Phenotype dataset ID: " << phenotype_dataset_id << endl;
  cout << "Phenotype file: " << phenotype_file << endl;
  cout << "Phenotype: " << phenotype << endl;
  for (auto&& kv : phenotype_column_types) {
    cout << ".. " << "column " << kv.first << " of type " << kv.second << endl;
  }
  cout << "Phenotype analysis columns: " << endl;
  for (auto&& v : phenotype_analysis_columns) {
    cout << ".. " << v << endl;
  }
  cout << "# rows: " << phenotype_nrows << endl;
  cout << "Sample subset key: " << sample_subset << endl;
  cout << "Samples provided? " << !samples.empty() << endl;
  cout << "Masks: " << endl;
  if (masks.empty()) {
    cout << ".. ** NO MASKS FOUND **" << endl;
  }
  else {
    for (auto&& m : masks) {
      m.print_groups(5, 5);
    }
  }

  cout << "Summary statistic dataset ID: " << summary_stat_dataset_id << endl;
  cout << "Score statistic files: " << endl;
  if (summary_stat_score_files.empty()) {
    cout << ".. ** NO SCORE STAT FILES FOUND ** " << endl;
  }
  else {
    for (auto&& f : summary_stat_score_files) {
      cout << ".. " << f << endl;
    }
  }

  cout << "Covariance files: " << endl;
  if (summary_stat_cov_files.empty()) {
    cout << ".. ** NO COV FILES FOUND ** " << endl;
  }
  else {
    for (auto&& f : summary_stat_cov_files) {
      cout << ".. " << f << endl;
    }
  }

  cout << "Segment size: " << segment_size << endl;
  cout << "Redis hostname: " << redis_hostname << endl;
  cout << "Redis port: " << redis_port << endl;
}

shared_ptr<ScoreCovarianceConfig> make_score_covariance_config() {
  return make_shared<ScoreCovarianceConfig>();
}

ScoreCovarianceRunner::ScoreCovarianceRunner(std::shared_ptr<ScoreCovarianceConfig> config) : config(config) {
  // Perform some basic checking of the configuration object to make sure we can actually complete this run.
  if (config->chrom.empty()) { throw std::invalid_argument("Must provide chromosome"); }
  if (config->start <= 0) { throw std::invalid_argument("Invalid starting position " + to_string(config->start)); }
  if (config->segment_size <= 0) { throw std::invalid_argument("Segment size must be non-zero"); }

  if (config->genotype_files.empty() && config->summary_stat_score_files.empty()) {
    throw std::invalid_argument("Must provide either genotype/phenotype files, or score stat/covariance files");
  }

  // If the config specifies a genotype file or phenotype file, we'll assume we're in "compute mode" - calculate scores
  // and covariances from the genotype & phenotype.
  if (!config->genotype_files.empty() || !config->phenotype_file.empty()) {
    run_mode = ScoreCovRunMode::COMPUTE;

    if (config->phenotype_file.empty()) {
      throw std::invalid_argument("Must provide phenotype file when genotype files are given");
    }
    else if (config->genotype_files.empty()) {
      throw std::invalid_argument("Must provide genotype files when a phenotype file is given");
    }
  }

  // If the config has a score statistic or covariance file, we're in "pre-compute mode" - the scores and covariances
  // are already computed, and we just need to read them out of files on disk.
  if (!config->summary_stat_score_files.empty() || !config->summary_stat_cov_files.empty()) {
    run_mode = ScoreCovRunMode::PRECOMPUTE;

    if (config->summary_stat_score_files.empty()) {
      throw std::invalid_argument("Must provide score statistic file in addition to covariance file");
    }
    else if (config->summary_stat_cov_files.empty()) {
      throw std::invalid_argument("Must provide covariance file in addition to score statistic file");
    }
  }

  if (run_mode == ScoreCovRunMode::COMPUTE) {
    // We're in LDServer/ScoreServer mode. In this mode, we compute score stats and covariance on-the-fly
    // from genotype and phenotype files.
    ld_server = make_shared<LDServer>(config->segment_size);
    score_server = make_shared<ScoreServer>(config->segment_size);

    for (auto&& genotype_file : config->genotype_files) {
      if (!is_file(genotype_file)) {
        throw std::invalid_argument("Genotype file not accessible: " + genotype_file);
      }

      ld_server->set_file(genotype_file);
      score_server->set_genotypes_file(genotype_file, config->genotype_dataset_id);
    }

    auto analysis_cols = make_shared<vector<string>>(config->phenotype_analysis_columns.begin(), config->phenotype_analysis_columns.end());

    score_server->load_phenotypes_file(
      config->phenotype_file,
      config->phenotype_column_types,
      config->phenotype_nrows,
      config->phenotype_delim,
      config->phenotype_sample_column,
      config->phenotype_dataset_id,
      analysis_cols
    );
    score_server->set_phenotype(config->phenotype);

    coordinate_samples(*score_server, *ld_server, config->genotype_files[0], config->phenotype, config->sample_subset, config->samples);

    //  if (!config->redis_hostname.empty()) {
    //    ld_server.enable_cache(config->genotype_dataset_id, config->redis_hostname, config->redis_port);
    //    score_server.enable_cache(config->redis_hostname, config->redis_port);
    //  }
  }
  else {
    // We have pre-computed score/covariance and just need to serve them from files.
    if (config->summary_stat_format == "METASTAAR") {
      summary_stat_loader = make_shared<MetastaarSummaryStatisticsLoader>(config->summary_stat_score_files, config->summary_stat_cov_files);
    }
    else {
      // Default format assumed to be raremetal/rvtest.
      summary_stat_loader = make_shared<RaremetalSummaryStatisticsLoader>(config->summary_stat_score_files, config->summary_stat_cov_files);
    }
  }
}

// TODO: fix copies below, make shared_ptr
void ScoreCovarianceRunner::run() {
  #ifndef NDEBUG
  cout << "Beginning run for configuration -- " << endl;
  config->pprint();
  #endif

  document = make_shared<Document>();
  document->SetObject();
  Document::AllocatorType& alloc = document->GetAllocator();

  Value data(kObjectType);
  Value variants(kArrayType);
  Value groups(kArrayType);

  set<string> seen_variants;
  auto ld_res = make_shared<LDQueryResult>(INITIAL_RESULT_SIZE);
  ld_res->limit = MAX_UINT32;

  auto score_res = make_shared<ScoreStatQueryResult>(INITIAL_RESULT_SIZE);
  score_res->limit = MAX_UINT32;

  for (auto&& mask : config->masks) {
    #ifndef NDEBUG
    cout << "Working on: " << mask.get_id() << endl;
    #endif
    for (auto&& group_item : mask) {
      #ifndef NDEBUG
      cout << ".. group: " << group_item.first << endl;
      #endif

      // Clear result objects
      ld_res->erase();
      score_res->erase();

      const VariantGroup& group = group_item.second;

      if (run_mode == ScoreCovRunMode::COMPUTE) {
        // In this block we use the LDServer and ScoreServer to compute scores and covariances from genotype/phenotype.
        SharedSegmentVector segments = make_shared_segment_vector();

        auto group_positions = group.get_positions();
        for_each(
          group_positions->begin(),
          group_positions->end(),
          [this](const uint64_t& p) {
            ld_server->add_overlap_position(p);
          }
        );

        ld_server->compute_region_ld(
          group.chrom,
          group.start,
          group.stop,
          correlation::COV,
          *ld_res,
          config->sample_subset,
          true, // compute diagonal elements (variance of each variant)
          segments
        );

        score_server->compute_scores(
          group.chrom,
          group.start,
          group.stop,
          *score_res,
          config->sample_subset,
          segments
        );
      }
      else {
        // Here we'll use the SummaryStatisticsLoader to read already computed scores/covariances from files on disk.
        try {
          summary_stat_loader->load_region(group.chrom, group.start, group.stop);
          ld_res = summary_stat_loader->getCovResult();
          score_res = summary_stat_loader->getScoreResult();
        }
        catch(NoVariantsInRange& e) {
          // There were no variants within the region requested, we can skip this group.
          continue;
        }
      }

      // A set of variants was explicitly asked for, so we filter down to them first.
      auto gvars = group.get_variants();
      if (!gvars->empty()) {
        ld_res->filter_by_variants(*group.get_variants());
        score_res->filter_by_variants(*group.get_variants());
      }

      // Run each user-defined filter passed in through the API. These can filter on MAF, p-value, or other properties
      // known to the score result object(s).
      for (auto& f : group.filters) {
        score_res->filter(f);
      }
      ld_res->filter_by_variants(*score_res->get_variants());

      if (score_res->data.size() == 0) {
        // No variants in this group survived after processing samples for missing data (they became monomorphic and
        // therefore could not be tested.
        continue;
      }

      score_res->sort_by_variant();

      for (auto&& v : score_res->data) {
        // The following code prevents duplicate variants from being inserted into the JSON.
        if (seen_variants.find(v.variant) == seen_variants.end()) {
          seen_variants.emplace(v.variant);
        }
        else {
          continue;
        }

        // TODO: internally at some point we should convert all variants stored in ScoreResult and related objects
        // to a VariantMeta object
        VariantMeta vm(v.variant);
        std::string variant_formatted;
        if (config->variant_format == VariantFormat::COLONS) {
          variant_formatted = vm.as_colons();
        }
        else if (config->variant_format == VariantFormat::EPACTS) {
          variant_formatted = vm.as_epacts();
        }
        else {
          throw std::invalid_argument("Bad variant format given when converting variant IDs");
        }

        Value this_variant(kObjectType);
        this_variant.AddMember("variant", Value(variant_formatted.c_str(), alloc), alloc);
        this_variant.AddMember("altFreq", v.alt_freq, alloc);

        std::isnan(v.pvalue) ? this_variant.AddMember("pvalue", Value(), alloc) : this_variant.AddMember("pvalue", v.pvalue, alloc);
        std::isnan(v.score_stat) ? this_variant.AddMember("score", Value(), alloc) : this_variant.AddMember("score", v.score_stat, alloc);
        variants.PushBack(this_variant, alloc);
      }

      Value this_group(kObjectType);
      this_group.AddMember("mask", mask.get_id(), alloc);
      this_group.AddMember("group", Value(group.name.c_str(), alloc), alloc);

      switch (mask.get_group_type()) {
        case VariantGroupType::GENE:
          this_group.AddMember("groupType", Value("GENE", alloc), alloc);
          break;
        case VariantGroupType::REGION:
          this_group.AddMember("groupType", Value("REGION", alloc), alloc);
      }

      Value group_variants(kArrayType);
      Value group_covar(kArrayType);
      set<string> seen_group_variants;

      for (unsigned int index1 = 0u; index1 < ld_res->data.variants.size(); ++index1) {
        auto& variant1 = ld_res->data.variants[index1];
        auto& offset1 = ld_res->data.offsets[index1];
        auto& correlations1 = ld_res->data.correlations[index1];

        for (unsigned int i = 0u; i < correlations1.size(); ++i) {
            auto& variant2 = ld_res->data.variants[i + offset1];
            double cov_value = correlations1[i];

            if (run_mode == ScoreCovRunMode::COMPUTE) {
                // When we compute from genotype/phenotype files, the covariance hasn't been normalized by the
                // residual variance yet, so we need to do so here. In the covariance matrix files, this has already been done.
                cov_value /= score_res->sigma2;
            }
            group_covar.PushBack(cov_value, alloc);

            if (seen_group_variants.find(variant1) == seen_group_variants.end()) {
                string variant1_formatted;
                if (config->variant_format == VariantFormat::COLONS) {
                    variant1_formatted = VariantMeta(variant1).as_colons();
                }
                else if (config->variant_format == VariantFormat::EPACTS) {
                    variant1_formatted = VariantMeta(variant1).as_epacts();
                }

                group_variants.PushBack(Value(variant1_formatted.c_str(), alloc), alloc);
                seen_group_variants.emplace(variant1);
            }

            if (seen_group_variants.find(variant2) == seen_group_variants.end()) {
                string variant2_formatted;
                if (config->variant_format == VariantFormat::COLONS) {
                    variant2_formatted = VariantMeta(variant2).as_colons();
                }
                else if (config->variant_format == VariantFormat::EPACTS) {
                    variant2_formatted = VariantMeta(variant2).as_epacts();
                }

                group_variants.PushBack(Value(variant2_formatted.c_str(), alloc), alloc);
                seen_group_variants.emplace(variant2);
            }
        }
      }

      this_group.AddMember("variants", group_variants, alloc);
      this_group.AddMember("covariance", group_covar, alloc);
      groups.PushBack(this_group, alloc);
    }
  }

  data.AddMember("variants", variants, alloc);
  data.AddMember("groups", groups, alloc);
  double& sigma2 = score_res->sigma2;
  double& nsamples = score_res->nsamples;
  std::isnan(sigma2) ? data.AddMember("sigmaSquared", Value(), alloc) : data.AddMember("sigmaSquared", sigma2, alloc);
  std::isnan(nsamples) ? data.AddMember("nSamples", Value(), alloc) : data.AddMember("nSamples", nsamples, alloc);

  if (run_mode == ScoreCovRunMode::COMPUTE) {
    data.AddMember("phenotypeDataset", config->phenotype_dataset_id, alloc);
    data.AddMember("genotypeDataset", config->genotype_dataset_id, alloc);
    data.AddMember("phenotype", Value(config->phenotype.c_str(), alloc), alloc);
  }
  else if (run_mode == ScoreCovRunMode::PRECOMPUTE) {
    data.AddMember("summaryStatisticDataset", config->summary_stat_dataset_id, alloc);
  }

  document->AddMember("data", data, alloc);
}

string ScoreCovarianceRunner::getJSON() const {
  StringBuffer strbuf;
  Writer<StringBuffer> writer(strbuf);

  // Note to future self: it is pretty common when testing on new data for Accept() to fail without a clear explanation
  // as to why. When that happens, breakpointing here and stepping into the function with the debugger helps figure out
  // what value the library is choking on. Frequently it is a NaN where one wasn't expected.
  if (!document->Accept(writer)) {
    throw runtime_error("Error while saving to JSON");
  }

  return strbuf.GetString();
}

string ScoreCovarianceRunner::getPrettyJSON() const {
  StringBuffer strbuf;
  PrettyWriter<StringBuffer> writer(strbuf);
  writer.SetIndent(' ', 2);
  if (!document->Accept(writer)) {
    throw runtime_error("Error while saving to JSON");
  }

  return strbuf.GetString();
}