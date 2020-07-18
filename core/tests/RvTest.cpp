#include "RvTest.h"
using namespace std;

void RvTestScores::load(const string &path) {
  unique_ptr<istream> file;
  ifstream fs(path, ios_base::in | ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;

  bool is_gz = path.find(".gz") != string::npos;
  if (is_gz) {
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(fs);
    file = make_unique<istream>(&inbuf);
  }
  else {
    file = make_unique<ifstream>(path);
  }

  string line;
  auto line_separator = regex("[ \t]");

  // Line regexes
  auto regex_samples = regex("##Samples=(\\d+)");
  auto regex_sigma = regex("## - Sigma2\t(.+)");
  auto regex_trait_sum = regex("##TraitSummary");
  auto regex_header = regex("CHROM\tPOS.*");

  bool header_done = false;
  bool parse_trait = false;
  while (getline(*file, line)) {
    smatch match;
    if (!header_done) {
      if (regex_search(line, match, regex_samples) && match.size() > 1) {
        this->nsamples = stoul(match.str(1));
      }
      else if (regex_search(line, match, regex_sigma) && match.size() > 1) {
        this->sigma2 = stod(match.str(1));
      }
      else if (regex_search(line, match, regex_trait_sum)) {
        parse_trait = true;
      }
      else if (parse_trait) {
        parse_trait = false;
        vector<string> trait_tok;
        copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(trait_tok));
      }
      else if (regex_search(line, match, regex_header)) {
        header_done = true;
      }
    }
    else {
      // Begin parsing record row
      vector<string> tokens;
      copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

      // For some reason RAREMETAL puts a genomic control line all the way at the end of the file
      // rvtest shouldn't do this, but it doesn't hurt to leave this check in anyway
      if (line.substr(0,1) == "#") { continue; }

      // Create record
      auto rec = make_shared<RvTestRecord>();
      rec->chrom = tokens.at(0);
      rec->pos = stoul(tokens.at(1));
      rec->ref = tokens.at(2);
      rec->alt = tokens.at(3);
      rec->n_informative = stoul(tokens.at(4));
      rec->founder_af = stod(tokens.at(5));
      rec->all_af = stod(tokens.at(5));
      rec->informative_alt_ac = stoul(tokens.at(6));
      rec->call_rate = stod(tokens.at(7));
      rec->hwe_pvalue = stod(tokens.at(8));
      rec->n_ref = stoul(tokens.at(9));
      rec->n_het = stoul(tokens.at(10));
      rec->n_alt = stoul(tokens.at(11));
      rec->u_stat = stod(tokens.at(12));
      rec->sqrt_vstat = stod(tokens.at(13));
      rec->alt_effsize = stod(tokens.at(14));
      rec->pvalue = stod(tokens.at(15));

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

RvTestScores::RvTestScores(const string &file) {
  load(file);
}

double RvTestScores::get_nsamples() {
  return nsamples;
}

double RvTestScores::get_sigma2() {
  return sigma2;
}

shared_ptr<RvTestRecord> RvTestScores::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<RvTestRecord>();
    return null;
  }
}