#ifndef LDSERVER_SCORESEGMENT_H
#define LDSERVER_SCORESEGMENT_H

#include "Segment.h"
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>

class ScoreSegment : public Segment {
protected:
  shared_ptr<vector<ScoreResult>> score_results;

public:
  ScoreSegment(const string& chromosome, uint64_t start_bp, uint64_t stop_bp, genotypes_store store);

  /**
   * Construct a new ScoreSegment, moving the data from an instance of the base class Segment.
   * @param other A Segment object.
   */
  ScoreSegment(Segment&& other) noexcept;
  bool has_scores() const;
  void compute_scores(const arma::vec& phenotype);
  void extract(uint64_t start, uint64_t end, struct ScoreStatQueryResult& result) const;
  void add_score(ScoreResult score);

  /**
   * Comparisons
   */
  bool operator==(const ScoreSegment& other) const;

  /**
   * Load/save functions for redis.
   * These are also overloaded below for loading/saving from serialized binary.
   *
   * TODO: these functions couldn't be used directly from the base class, because the load/save functions
   * that take archive parameters can't be virtualed (because they are templated.) There's probably a better way
   * to do this but for now it's just directly copied code.
   *
   * @param redis_cache
   * @param key
   */
  virtual void load(redisContext* redis_cache, const string& key) override;
  virtual void save(redisContext* redis_cache, const string& key) override;

  /**
   * Load/save functions for binary format.
   * Stores the same elements as the base class Segment, and additionally score stats/pvalues/etc.
   *
   * Cereal's docs seem to think it's fine to hide the non-virtual method of the base class, since it can't be
   * virtual in the first place (can't have templated virtual functions.) If there's a better way to do this,
   * it would be nice.
   *
   * @tparam Archive
   * @param ar
   */
  template <class Archive> void load(Archive& ar) {
    ar(cereal::base_class<Segment>(this), score_results);
  }

  template <class Archive> void save(Archive& ar) const {
    ar(cereal::base_class<Segment>(this), score_results);
  }
};

#endif //LDSERVER_SCORESEGMENT_H
