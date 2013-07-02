#include <VecbosEGObject.hh>
#include <HggEGEnergyCorrector.hh>
#include <HggVertexing.hh>
#include <HggEnergyScale.hh>

#include <vector>
#include <iostream>
#include <string>

#include <TChain.h>
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TH1F.h"

using namespace std;
#include "HggVertexing.hh"

class ZeeSelector{
public:
  ZeeSelector();
  ~ZeeSelector();
  ZeeSelector(vector<string> fNames,string treeName,string outputFile);
  void loadChain(vector<string> fNames, string treeName);
  void setOutputFile(string s){outputFile = s;}
  bool isValid(){return valid;}
  void setIsData(bool d){isData_=d;}
  void setConfig(string s){config=s;}
  void Loop();
private:
  bool valid;;
  TChain* fChain;
  TTree* outTree;
  TTree* outTreeZee;
  string outputFile;
  string config;

  bool isData_;

  int init();
  void setBranchAddresses();
  void setupOutputTree();
  
  bool passPresel(VecbosEle&);

  HggEGEnergyCorrector *elecorr;

  // Mass Selection
  float DZmassref; 
  float DZmass;
  float Zeemass;
  int   nEle;
  float lpass;
  float tpass;
  float mvapass;
  bool  isZmass;
  float rho;
  float PFIsoOverPT1;
  float PFIsoOverPT2;

  // Input variables
  int runNumber;
  int evtNumber;
  int nVtx;

  // Electron Selection
  std::vector<VecbosEle> *Electrons;

  // Variables that will be outputted
  float mass;
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
