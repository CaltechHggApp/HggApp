#ifndef VecbosJetCorrector_hh
#define VecbosJetCorrector_hh

#include "FactorizedJetCorrector.h"
#include "JetCorrectorParameters.h"
#include "VecbosEGObject.hh"
#include "ReadConfig.hh"
#include <string>

class VecbosJetCorrector{
public:
  VecbosJetCorrector();
  VecbosJetCorrector(std::string);
  VecbosJetCorrector(ReadConfig&);
  void getConfig(std::string);
  void getConfig(ReadConfig&);
  void CorrectJet(VecbosJet&,float);
private:
  std::vector<JetCorrectorParameters> vPar;
  FactorizedJetCorrector *jfCorr;
};

#endif
