#ifndef LDSERVER_RAREMETALRUNNER_H
#define LDSERVER_RAREMETALRUNNER_H

#include <memory>
#include <string>
#include <vector>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/prettywriter.h>
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
  std::shared_ptr<rapidjson::Document> document;
public:
  void operator()(const std::vector<Mask>& masks, const std::string& sample_subset, const ScoreServer& score_server, const LDServer& ld_server);
  string getJSON() const;
  string getPrettyJSON() const;
};

#endif //LDSERVER_RAREMETALRUNNER_H
