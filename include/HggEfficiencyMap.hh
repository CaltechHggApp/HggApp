#include <VecbosEGObject.hh>
#include <HggPhotonID.hh>
#include <HggMuonID.hh>

#include <HggVertexing.hh>

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

class HggEfficiencyMap{
public:
  HggEfficiencyMap();
  ~HggEfficiencyMap();
  HggEfficiencyMap(vector<string> fNames,string treeName,string outputFile);
  void loadChain(vector<string> fNames, string treeName);
  void setOutputFile(string s){outputFile = s;}
  void setConfigFile(string s){configFile = s;}
  void setMuMuGamma(){mode=kMMG;}
  void setMC(){mode=kMC;}
  void setZeeMC(){mode=kZeeMC;}
  void setZeeData(){mode=kZee;}

  void Loop();
protected:
  TChain* fChain;
  bool valid;

  string outputFile;

  //Photon ID
  string configFile;
  HggPhotonID *PhotonID;

  //selection
  float leadPhoEtMin;
  float subleadPhoEtMin;
  bool doPtOverM;
  float PtOverMLead;
  float PtOverMSubLead;
  const static float rhoFac    = 0.17;
  const static float rhoFacBad = 0.52;
  const static float isoSumConst = 0;//5;
  const static float isoSumConstBad = 0;//7;

  int init();
  void clear();
  void setBranchAddresses();
  void setupOutputTree();

  bool etaSelectPhoton(VecbosPho* pho);
  void smearPhoton(VecbosPho* pho, int smearShift);

  enum runMode{kMC,kMMG,kZeeMC,kZee};
  runMode mode; //!< This flag will identify photons using MuMuGamma Selection in data

  void runGen();
  void runMMG();
  void runZeeMC();
  void runZee();

  VecbosPho* findGenMatchPhoton(VecbosGen* gen); //find the reco photon matched to a gen photon
  VecbosPho* findGenMatchElectron(VecbosGen* gen); //find the reco photon matched to a gen photon

  //*****************
  //output tree
  //****************
  TTree *outputTree;
  void setOutputBranches();

  //output variables:
  float outPt;
  float outRecoPt;
  float outEta;
  float outRecoEta;
  float outEtaSC;
  float outRecoEtaSC;
  float outPhi;
  float outR9;
  bool passMatch;
  bool passEta;
  bool passPreselection;
  bool passCiC;

  //***********
  //member data
  //***********
  const static int maxPho = 100;
  int nPho_;
  std::vector<VecbosPho> *Photons_; // this contains ALL photons
  
  bool photonMatchedElectron[maxPho];

  int nMu_;
  std::vector<VecbosMu> *Muons_;

  int nJet_;
  std::vector<VecbosJet> *Jets_;
  
  // for each photon, a vector of floats giving the track iso from each ggVertex

  int nVtx; 
  static const int MAXVX = 100;
  float vtxX[MAXVX];
  float vtxY[MAXVX];
  float vtxZ[MAXVX];
  float vtxChi2[MAXVX];
  float vtxNdof[MAXVX];
  float vtxNormalizedChi2[MAXVX];
  int vtxTrackSize[MAXVX];
  int vtxIsFake[MAXVX];
  int vtxIsValid[MAXVX];

  float rho;
  float rhoEtaMax44;

  int lumiBlock;
  int runNumber;
  long evtNumber;
  bool _isData;

  float evtWeight;
  float pileupWeight;
  
  int nGenHiggs;
  std::vector<VecbosGen> *GenHiggs;

  int nGenPho;
  std::vector<VecbosGen> *GenPhotons;

  int nGenEle;
  std::vector<VecbosGen> *GenElectrons;
  

  float pfMet;
  float pfMetPhi;

  float inPU;

  bool eeBadScFilterFlagOut;
  bool hcalLaserEventFilterFlagOut;
  bool HBHENoiseFilterResultFlagOut;
  bool isNotDeadEcalClusterOut;
  bool trackerFailureFilterFlagOut;
  bool CSCHaloFilterFlagOut;
  bool drDeadOut; 
  bool drBoundaryOut;
  bool ECALTPFilterFlagOut;

  bool eeBadScFilterFlag;
  bool hcalLaserEventFilterFlag;
  bool HBHENoiseFilterResultFlag;
  bool isNotDeadEcalCluster;
  bool trackerFailureFilterFlag;
  bool CSCHaloFilterFlag;
  bool drDead; 
  bool drBoundary;
  bool ECALTPFilterFlag;
};
