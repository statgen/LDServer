#ifndef LDSERVER_RAREMETAL_H
#define LDSERVER_RAREMETAL_H

#include <string>
#include <map>
#include <vector>
#include <memory>

using namespace std;

struct RareMetalRecord {
  string chrom;
  uint64_t pos;
  string ref;
  string alt;
  uint64_t n_informative;
  double founder_af;
  double all_af;
  uint64_t informative_alt_ac;
  double call_rate;
  double hwe_pvalue;
  uint64_t n_ref;
  uint64_t n_het;
  uint64_t n_alt;
  double u_stat;
  double sqrt_vstat;
  double alt_effsize;
  double pvalue;
};

class RareMetalScores {
protected:
  uint64_t nsamples;
  double sigma2; // sigma_e2_hat from RAREMETAL
  string trait_name;
  vector<shared_ptr<RareMetalRecord>> records;
  map<string, shared_ptr<RareMetalRecord>> index;
public:
  RareMetalScores(const string &file);
  void load(const string &file);
  double get_sigma2();
  double get_nsamples();
  shared_ptr<RareMetalRecord> get_record(const string &i);
};

#endif //LDSERVER_RAREMETAL_H
