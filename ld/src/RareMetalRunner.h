#ifndef LDSERVER_RAREMETALRUNNER_H
#define LDSERVER_RAREMETALRUNNER_H

#include <string>
#include <vector>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/writer.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include "Mask.h"
#include "LDServer.h"
#include "ScoreServer.h"
#include "Phenotypes.h"
#include "Segment.h"
#include "Types.h"

class RareMetalRunner {
protected:
  std::string json;
public:
  void operator()(const std::vector<Mask>& masks, const std::string& sample_subset, const ScoreServer& score_server, const LDServer& ld_server);
  string getJSON() const;
};

#endif //LDSERVER_RAREMETALRUNNER_H
