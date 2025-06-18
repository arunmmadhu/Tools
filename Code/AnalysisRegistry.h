#ifndef ANALYSISREGISTRY_H
#define ANALYSISREGISTRY_H

#include <map>
#include <string>
#include "Selection_Base.h"  // for Selection_Base*

class AnalysisRegistry {
public:
  static AnalysisRegistry& instance() {
    static AnalysisRegistry inst;
    return inst;
  }

  void registerAnalysis(const std::string& name, Selection_Base* a) {
    map_[name] = a;
  }

  Selection_Base* get(const std::string& name) const {
    auto it = map_.find(name);
    return it != map_.end() ? it->second : nullptr;
  }

private:
  std::map<std::string, Selection_Base*> map_;
};

#endif
