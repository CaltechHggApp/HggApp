#ifndef SusyHggSelectorSimple_h
#define SusyHggSelectorSimple_h

#include "SusyHggSelector.hh"
#include "BaseSelector.hh"
#include "StandardPhotonID.hh"
#include "StandardElectronID.hh"
#include "VecbosJetID.hh"
#include "RazorVariables.hh"
#include "JECUReader.hh"
#include "HggEnergyScale.hh"


#include "TVector3.h"
#include "TLorentzVector.h"

#include "assert.h"
#include <bitset>

class SusyHggSelectorSimple : public SusyHggSelector {
public:
  SusyHggSelectorSimple(std::vector<std::string> fNames, std::string treeName,std::string outputFile);
  void setScalingCFG(std::string scaleCFG);
protected:

  HggEnergyScale  *scaler=0;
  //mandatory overrides
  virtual void processEntry(Long64_t iEntry);

  virtual void clear();
};

#endif
