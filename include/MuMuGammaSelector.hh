#include <VecbosEGObject.hh>
#include <HggEGEnergyCorrector.hh>
#include <HggEnergyScale.hh>
#include <HggPhotonID.hh>

#include <vector>
#include <iostream>
#include <string>

#include <TChain.h>
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TH1F.h"

#include "BaseSelector.hh"

using namespace std;
//#include "HggVertexing.hh"

class MuMuGammaSelector : public BaseSelector{
public:
  MuMuGammaSelector(vector<string> fNames,string treeName,string outputFile):BaseSelector(fNames,treeName,outputFile){}
protected:
  HggPhotonID *photonID;

  //mandatory overrides
  virtual void processEntry(Long64_t iEntry);
  virtual int init(){ return 0; }
  virtual void processConfig(ReadConfig &cfg);
  virtual void setupOutputTree();
  
  virtual void clear();
  virtual void firstInit();

  bool passPresel(VecbosMu&);
  void writeEventInfo();


  //output variables
  int runNumberOut;
  int evtNumberOut;
  int nVtxOut;
  int rhoOut;

  const static int maxMuMuG=200;
  int nMuMuG;
  float massMuMuGamma[maxMuMuG];
  float massMuMuRegGamma[maxMuMuG];
  float massMuMuScaleGamma[maxMuMuG];
  float massMuMuGenGamma[maxMuMuG];
  float massMuMu[maxMuMuG];
  float puWeight[maxMuMuG];
  
  std::vector<VecbosMu> outMuon1;
  std::vector<VecbosMu> outMuon2;
  std::vector<VecbosPho> outPhoton;
  float isosumoetPho[maxMuMuG];
  float mvaPho[maxMuMuG];

  
};
