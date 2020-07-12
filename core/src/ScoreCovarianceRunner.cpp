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

  if (config->genotype_files.empty() && config->summary_stat_score_file.empty()) {
    throw std::invalid_argument("Must provide either genotype/phenotype files, or score stat/covariance files");
  }

  if (!config->genotype_files.empty() || !config->phenotype_file.empty()) {
    run_mode = ScoreCovRunMode::COMPUTE;

    if (config->phenotype_file.empty()) {
      throw std::invalid_argument("Must provide phenotype file when genotype files are given");
    }
    else if (config->genotype_files.empty()) {
      throw std::invalid_argument("Must provide genotype files when a phenotype file is given");
    }
  }

  if (!config->summary_stat_score_file.empty() || !config->summary_stat_cov_file.empty()) {
    run_mode = ScoreCovRunMode::PRECOMPUTE;

    if (config->summary_stat_score_file.empty()) {
      throw std::invalid_argument("Must provide score statistic file in addition to covariance file");
    }
    else if (config->summary_stat_cov_file.empty()) {
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
    summary_stat_loader = make_shared<SummaryStatisticsLoader>(config->summary_stat_score_file, config->summary_stat_cov_file);
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
        summary_stat_loader->load_region(group.chrom, group.start, group.stop);
        ld_res = summary_stat_loader->getCovResult();
        score_res = summary_stat_loader->getScoreResult();
      }

      ld_res->filter_by_variants(*group.get_variants());
      score_res->filter_by_variants(*group.get_variants());

      if (score_res->data.size() == 0) {
        // No variants in this group survived after processing samples for missing data (they became monomorphic and
        // therefore could not be tested.
        continue;
      }

      ld_res->sort_by_variant();
      score_res->sort_by_variant();

      for (auto&& v : score_res->data) {
        if (seen_variants.find(v.variant) == seen_variants.end()) {
          seen_variants.emplace(v.variant);
        }
        else {
          continue;
        }

        Value this_variant(kObjectType);
        this_variant.AddMember("variant", Value(v.variant.c_str(), alloc), alloc);
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
      for (auto&& pair : ld_res->data) {
        double cov_value = pair.value;
        if (run_mode == ScoreCovRunMode::COMPUTE) {
          // When we compute from genotype/phenotype files, the covariance hasn't been normalized by the
          // residual variance yet, so we need to do so here. In the covariance matrix files, this has already been done.
          cov_value /= score_res->sigma2;
        }
        group_covar.PushBack(cov_value, alloc);

        if (seen_group_variants.find(pair.variant1) == seen_group_variants.end()) {
          group_variants.PushBack(Value(pair.variant1.c_str(), alloc), alloc);
          seen_group_variants.emplace(pair.variant1);
        }

        if (seen_group_variants.find(pair.variant2) == seen_group_variants.end()) {
          group_variants.PushBack(Value(pair.variant2.c_str(), alloc), alloc);
          seen_group_variants.emplace(pair.variant2);
        }
      }

      this_group.AddMember("variants", group_variants, alloc);
      this_group.AddMember("covariance", group_covar, alloc);
      groups.PushBack(this_group, alloc);
    }
  }

  data.AddMember("variants", variants, alloc);
  data.AddMember("groups", groups, alloc);
  data.AddMember("sigmaSquared", score_res->sigma2, alloc);
  data.AddMember("nSamples", score_res->nsamples, alloc);
  data.AddMember("phenotypeDataset", config->phenotype_dataset_id, alloc);
  data.AddMember("genotypeDataset", config->genotype_dataset_id, alloc);
  data.AddMember("phenotype", Value(config->phenotype.c_str(), alloc), alloc);
  document->AddMember("data", data, alloc);
}

string ScoreCovarianceRunner::getJSON() const {
  StringBuffer strbuf;
  Writer<StringBuffer> writer(strbuf);
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