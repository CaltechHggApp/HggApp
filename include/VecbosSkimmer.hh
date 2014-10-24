#ifndef VecbosSkimmer_hh
#define VecbosSkimmer_hh

#include "VecbosBase.hh"
#include "PassTrigger.hh"

#include <unordered_set>
#include "stdint.h"

class VecbosSkimmer : public VecbosBase {
public:
  VecbosSkimmer(TTree* tree) : VecbosBase(tree) {}

  void addEvent(uint32_t run, uint32_t lumi, uint32_t event) {
    eventsToSkim.insert(evtInfo(run,lumi,event));
  }
  void Skim(const char* outputFileName);

private:
  std::unordered_set<evtInfo, eiHash> eventsToSkim;
};

#endif
