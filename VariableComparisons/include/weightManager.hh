#ifndef weightManager_hh
#define weightManager_hh

#include "TString.h"
#include <map>

class weightManager{
public:
  weightManager();
  float getWeight(TString name,float lumi=1); //lumi in 1/fb
protected:
  std::map<TString,float> Nevents;
  std::map<TString,float> crossSection;
};

#endif
