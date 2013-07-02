#include <VecbosEGObject.hh>
#include <HggEGEnergyCorrector.hh>
#include <HggVertexing.hh>
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

using namespace std;
#include "HggVertexing.hh"

class MuMuGammaSelector{
public:
  //MuMuGammaSelector();
  ~MuMuGammaSelector();
  MuMuGammaSelector(vector<string> fNames,string treeName,string outputFile);
  void loadChain(vector<string> fNames, string treeName);
  void setOutputFile(string s){outputFile = s;}
  bool isValid(){return valid;}
  void setIsData(bool d){isData_=d;}
  void Loop();
  void setConfigFile(string s){cfg=s;}
private:
  bool valid;;
  TChain* fChain;
  TTree* outTree;
  string outputFile;

  string cfg;

  bool isData_;

  HggPhotonID *photonID;

  int init();
  void setBranchAddresses();
  void setupOutputTree();
  
  void clear();
  bool passPresel(VecbosMu&);
  void writeEventInfo();

  // Input variables
  int runNumber;
  int evtNumber;
  int rho;

  //vtx
  static const int maxVtx=100;
  
  int nVtx;
  float vtxX[maxVtx];
  float vtxY[maxVtx];
  float vtxZ[maxVtx];

  // Electron Selection
  int __nMu;
  MuCollection  *__Muons;
  int __nPho;
  PhoCollection *__Photons;

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
