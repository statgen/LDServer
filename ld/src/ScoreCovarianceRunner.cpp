#include "ScoreCovarianceRunner.h"
using namespace std;
using namespace rapidjson;

const uint32_t MAX_UINT32 = numeric_limits<uint32_t>::max();
const uint32_t INITIAL_RESULT_SIZE = 100000000;

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
  cout << "# rows: " << phenotype_nrows << endl;
  cout << "Sample subset key: " << sample_subset << endl;
  cout << "Samples provided? " << !samples.empty() << endl;
  cout << "Masks: " << endl;
  if (masks.empty()) {
    cout << ".. ** NO MASKS FOUND **" << endl;
  }
  else {
    for (auto&& m : masks) {
      cout << ".. " << m.get_id() << endl;
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
  if (config->phenotype_file.empty()) { throw std::invalid_argument("Must provide phenotype file"); }
  if (config->genotype_files.empty()) { throw std::invalid_argument("Must provide genotype file"); }
  if (config->segment_size <= 0) { throw std::invalid_argument("Segment size must be non-zero"); }
}

// TODO: fix copies below, make shared_ptr
void ScoreCovarianceRunner::run() {
  #ifndef NDEBUG
  cout << "Beginning run for configuration -- " << endl;
  config->pprint();
  #endif

  LDServer ld_server(config->segment_size);
  ScoreServer score_server(config->segment_size);

  for (auto&& genotype_file : config->genotype_files) {
    ld_server.set_file(genotype_file);
    score_server.set_genotypes_file(genotype_file, config->genotype_dataset_id);
  }

  if (config->sample_subset != "ALL") {
    ld_server.set_samples(config->sample_subset, config->samples);
    score_server.set_samples(config->sample_subset, config->samples);
  }

  score_server.load_phenotypes_file(
    config->phenotype_file,
    config->phenotype_column_types,
    config->phenotype_nrows,
    config->phenotype_delim,
    config->phenotype_sample_column,
    config->phenotype_dataset_id
  );
  score_server.set_phenotype(config->phenotype);

//  if (!config->redis_hostname.empty()) {
//    ld_server.enable_cache(config->genotype_dataset_id, config->redis_hostname, config->redis_port);
//    score_server.enable_cache(config->redis_hostname, config->redis_port);
//  }

  document = make_shared<Document>();
  document->SetObject();
  Document::AllocatorType& alloc = document->GetAllocator();

  Value data(kObjectType);
  Value variants(kArrayType);
  Value groups(kArrayType);

  set<string> seen_variants;
  LDQueryResult ld_res(INITIAL_RESULT_SIZE);
  ld_res.limit = MAX_UINT32;
  ScoreStatQueryResult score_res(INITIAL_RESULT_SIZE);
  score_res.limit = MAX_UINT32;
  for (auto&& mask : config->masks) {
    #ifndef NDEBUG
    cout << "Working on: " << mask.get_id() << endl;
    #endif
    for (auto&& group_item : mask) {
      #ifndef NDEBUG
      cout << ".. group: " << group_item.first << endl;
      #endif

      // Clear result objects
      ld_res.erase();
      score_res.erase();

      const VariantGroup& group = group_item.second;

      // DEBUG PATCH

      // END PATCH

      SharedSegmentVector segments = make_shared_segment_vector();

      ld_server.compute_region_ld(
        group.chrom,
        group.start,
        group.stop,
        correlation::COV,
        ld_res,
        config->sample_subset,
        true, // compute diagonal elements (variance of each variant)
        segments
      );

      score_server.compute_scores(
        group.chrom,
        group.start,
        group.stop,
        score_res,
        config->sample_subset,
        segments
      );

      ld_res.filter_by_variants(*group.get_variants());
      score_res.filter_by_variants(*group.get_variants());

      ld_res.sort_by_variant();
      score_res.sort_by_variant();

      for (auto&& v : score_res.data) {
        if (seen_variants.find(v.variant) == seen_variants.end()) {
          seen_variants.emplace(v.variant);
        }
        else {
          continue;
        }

        Value this_variant(kObjectType);
        this_variant.AddMember("variant", Value(v.variant.c_str(), alloc), alloc);
        this_variant.AddMember("altFreq", v.alt_freq, alloc);
        this_variant.AddMember("pvalue", v.pvalue, alloc);
        this_variant.AddMember("score", v.score_stat, alloc);
        variants.PushBack(this_variant, alloc);
      }

      Value this_group(kObjectType);
      this_group.AddMember("mask", Value(mask.get_id().c_str(), alloc), alloc);
      this_group.AddMember("group", Value(group.name.c_str(), alloc), alloc);

      switch (mask.get_group_type()) {
        case VariantGroupType::GENE:
          this_group.AddMember("groupType", Value("gene", alloc), alloc);
          break;
        case VariantGroupType::REGION:
          this_group.AddMember("groupType", Value("region", alloc), alloc);
      }

      Value group_variants(kArrayType);
      Value group_covar(kArrayType);
      set<string> seen_group_variants;
      for (auto&& pair : ld_res.data) {
        group_covar.PushBack(pair.value / score_res.sigma2, alloc);

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
  data.AddMember("sigmaSquared", score_res.sigma2, alloc);
  data.AddMember("nSamples", score_res.nsamples, alloc);
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