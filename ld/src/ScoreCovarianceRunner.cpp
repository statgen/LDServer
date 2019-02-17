#include "RareMetalRunner.h"
using namespace std;
using namespace rapidjson;

const uint32_t MAX_UINT32 = numeric_limits<uint32_t>::max();

// TODO: fix copies below, make shared_ptr
void RareMetalRunner::operator()(const vector<Mask>& masks, const string& sample_subset, const ScoreServer& score_server, const LDServer& ld_server) {
  document = make_shared<Document>();
  document->SetObject();
  Document::AllocatorType& alloc = document->GetAllocator();

  Value data(kObjectType);
  Value variants(kArrayType);
  Value groups(kArrayType);

  set<string> seen_variants;
  LDQueryResult ld_res(MAX_UINT32);
  ScoreStatQueryResult score_res(MAX_UINT32);
  for (auto&& mask : masks) {
    for (auto&& group_item : mask) {
      // Clear result objects
      ld_res.erase();
      score_res.erase();

      const VariantGroup& group = group_item.second;

      SharedSegmentVector segments = make_shared_segment_vector();

      ld_server.compute_region_ld(
        group.chrom,
        group.start,
        group.stop,
        correlation::COV,
        ld_res,
        sample_subset,
        segments
      );

      score_server.compute_scores(
        group.chrom,
        group.start,
        group.stop,
        score_res,
        sample_subset,
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
      this_group.AddMember("mask", Value(mask.get_name().c_str(), alloc), alloc);

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
        group_covar.PushBack(pair.value, alloc);

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

string RareMetalRunner::getJSON() const {
  StringBuffer strbuf;
  Writer<StringBuffer> writer(strbuf);
  if (!document->Accept(writer)) {
    throw runtime_error("Error while saving to JSON");
  }

  return strbuf.GetString();
}

string RareMetalRunner::getPrettyJSON() const {
  StringBuffer strbuf;
  PrettyWriter<StringBuffer> writer(strbuf);
  if (!document->Accept(writer)) {
    throw runtime_error("Error while saving to JSON");
  }

  return strbuf.GetString();
}