#include <VecbosEGObject.hh>
#include <HggEGEnergyCorrector.hh>
#include <HggEnergyScale.hh>

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

struct ElectronAdditionalInfo{
  float dEoE,dEoEErr;
  float scaledEnergy,scaledEnergyError;
};

class ZeeSelector : public BaseSelector{
public:
  ZeeSelector(vector<string> fNames,string treeName,string outputFile):BaseSelector(fNames,treeName,outputFile) {}
private:
  bool doSmear=false;
  bool doScale=false;
  int applyScaleSmear=0;

  HggEnergyScale *scale=0;
  HggEnergyScale *smear=0;

  virtual int init();
  virtual void setupOutputTree();
  virtual void processEntry(Long64_t iEntry);


  virtual void clear();
  virtual void firstInit();
  virtual void processConfig(ReadConfig &cfg);

  bool passPresel(VecbosEle&);



  HggEGEnergyCorrector *elecorr;
  std::vector<ElectronAdditionalInfo> eleInfo;


  // Mass Selection
  float DZmassref; 
  float Zeemass;
  float lpass;
  float tpass;
  float mvapass;
  bool  isZmass;
  float rho;
  float PFIsoOverPT1;
  float PFIsoOverPT2;

  // Variables that will be outputted
  float mass;
  float DZmass;
  int   nEleOut;
  float Ele1mva;
  float Ele2mva;

  float Ele1pt;
  float Ele1eta;
  float Ele1phi;
  float Ele1E;
  float Ele1Epho;

  float Ele2pt;
  float Ele2eta;
  float Ele2phi;
  float Ele2E;
  float Ele2Epho;

  float Ele1etaSC;
  float Ele2etaSC;

  float Ele1r9;
  float Ele2r9;
  float Ele1sigEoE;
  float Ele2sigEoE;
  float Ele1sigEoEpho;
  float Ele2sigEoEpho;
  float Ele1sigEscaleoEpho;
  float Ele2sigEscaleoEpho;
  float passtight;
  float passmva;	
  float passloose;
};
