#include "RareMetal.h"
#include <fstream>
#include <string>
#include <iostream>
#include <regex>

using namespace std;

void RareMetalScores::load(const string &file) {
  ifstream input_file(file);
  string line;
  auto line_separator = regex("[ \t]");

  // Line regexes
  auto regex_samples = regex("##Samples=(\\d+)");
  auto regex_sigma = regex("##Sigma_e2_Hat\t(.+)");
  auto regex_trait_sum = regex("##TraitSummaries");
  auto regex_header = regex("#CHROM\tPOS.*");

  bool header_done = false;
  bool parse_trait = false;
  while (getline(input_file, line)) {
    smatch match;
    if (!header_done) {
      if (regex_search(line, match, regex_samples) && match.size() > 1) {
        this->nsamples = stoul(match.str(1));
      }
      else if (regex_search(line, match, regex_sigma) && match.size() > 1) {
        this->sigma = stod(match.str(1));
      }
      else if (regex_search(line, match, regex_trait_sum)) {
        parse_trait = true;
      }
      else if (parse_trait) {
        parse_trait = false;
        vector<string> trait_tok;
        copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(trait_tok));
        this->trait_name = trait_tok.at(0);
      }
      else if (regex_search(line, match, regex_header)) {
        header_done = true;
      }
    }
    else {
      // Begin parsing record row
      vector<string> tokens;
      copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

      // For some reason RAREMETAL puts a genomic control line all the way at the end of the file...
      if (line.substr(0,1) == "#") { continue; }

      // Create record
      auto rec = make_shared<RareMetalRecord>();
      rec->chrom = tokens.at(0);
      rec->pos = stoul(tokens.at(1));
      rec->ref = tokens.at(2);
      rec->alt = tokens.at(3);
      rec->n_informative = stoul(tokens.at(4));
      rec->founder_af = stod(tokens.at(5));
      rec->all_af = stod(tokens.at(6));
      rec->informative_alt_ac = stoul(tokens.at(7));
      rec->call_rate = stod(tokens.at(8));
      rec->hwe_pvalue = stod(tokens.at(9));
      rec->n_ref = stoul(tokens.at(10));
      rec->n_het = stoul(tokens.at(11));
      rec->n_alt = stoul(tokens.at(12));
      rec->u_stat = stoul(tokens.at(13));
      rec->sqrt_vstat = stod(tokens.at(14));
      rec->alt_effsize = stod(tokens.at(15));
      rec->pvalue = stod(tokens.at(16));

      // Keys
      string chrpos = rec->chrom + ":" + to_string(rec->pos);
      string variant = rec->chrom + ":" + to_string(rec->pos) + "_" + rec->ref + "/" + rec->alt;

      // Insert
      records.push_back(rec);
      this->index.emplace(chrpos, rec);
      this->index.emplace(variant, rec);
    }
  }
}

RareMetalScores::RareMetalScores(const string &file) {
  load(file);
}

double RareMetalScores::get_nsamples() {
  return nsamples;
}

double RareMetalScores::get_sigma() {
  return sigma;
}

shared_ptr<RareMetalRecord> RareMetalScores::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<RareMetalRecord>();
    return null;
  }
}