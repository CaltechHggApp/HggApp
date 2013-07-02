#include <VecbosEGObject.hh>
#include <HggMassResolution.hh>
#include <HggPhotonID.hh>

#include <HggEGEnergyCorrector.hh>
//#include <HggVertexing.hh>
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
//#include "HggVertexing.hh"

class HggSelector{
public:
  HggSelector();
  ~HggSelector();
  HggSelector(vector<string> fNames,string treeName,string outputFile);
  void loadChain(vector<string> fNames, string treeName);
  void setOutputFile(string s){outputFile = s;}
  void setConfig(string s){configFile = s;}
  void setMassResConfig(string s){massResConfig = s;}
  bool isValid(){return valid;}
  void addTrigger(string s){triggers.push_back(s);}
  void setDoMuMuGamma(bool d=true){doMuMuGamma=d;}
  void setNSigma(int ns){nSigma=ns;}

  void suppressElectronVeto(){doElectronVeto=false;}
  void setForceVertexZero(){forceVtxZero=true;}
  void setIsData(bool d){isData_=d;}
  void Loop();
private:
  TChain* fChain;
  bool valid;
  bool doElectronVeto;
  bool doMuMuGamma;
  bool forceVtxZero;
  bool isData_;

  TTree* outTree;
  TTree* outTreeMuMuG;
  string configFile;
  string outputFile;

  std::pair<int,int> getBestPair(float*,int,int);
  std::pair<int,int> getBestPairCiC(int,int,bool);
  

  bool preSelectPhotons(VecbosPho*,VecbosPho*,TVector3); // kinematic photons selections
  
  string massResConfig;
  HggMassResolution *massRes;

  HggPhotonID *PhotonID;

  //incase we want to redo the regression/scaling/smearing
  bool doRegression;
  bool doScale;
  bool doSmear;
  HggEGEnergyCorrector* corrector;
  HggEnergyScale* scale;
  int applyScaleSmear;
  HggEnergyScale* smear;

  void smearPhoton(VecbosPho*,int);

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

  vector<string> triggers;
  int *triggerDec;
  bool requireTrigger();

  int getVertexIndex(int,int);
  float getVertexProb(int,int);
  float getVertexMVA(int,int);

  int init();
  void clear();
  void setDefaults();
  void setBranchAddresses();
  void setupOutputTree();
  void setupTMVA();
  void fillGenInfo();

  void fillMuMuGamma();
  float getVBFMjj(VecbosPho*,VecbosPho*,TVector3,float*);
  bool passJetID(VecbosJet*);
  //TMVA stuff
  string weightFile_diPho;
  string methodName_diPho;
  TMVA::Reader *diPhotonMVA;

  std::map<std::string,TH1F*> MVAInputs;
  std::map<std::string,TH1F*> MjjDists;
  float getDiPhoMVA(int,int,float,float,bool);

  float getMPair(int,int);

  ReducedPhotonData getReducedData(VecbosPho*,TVector3,int);
  //These are the variables that will be filled for the TMVA:
  //better to do it this way so we can do PFSC or regular SC without major changes

  
  //diPhotonMVA:
  float smearedMassErrByMass; // mass error smeared/mass
  float smearedMassErrByMassWrongVtx;
  float vtxprob;
  float pho1PtByMass;
  float pho2PtByMass;
  float pho1Eta;
  float pho2Eta;
  float cosDPhi;
  float pho1IdMVA;
  float pho2IdMVA;
  //-------------

  int getCategory();
  int getCategoryPFCiC();

  //Output Variables:
  int trigger_;
  float mPair_;
  float mPairNoCorr_;
  float mPairRes_;
  float mPairResWrongVtx_;
  float diPhoMVA_;
  float diPhoMVAShift_[9];
  int diPhoVtx_;
  float vtxProb_;
  float diPhoVtxX_;
  float diPhoVtxY_;
  float diPhoVtxZ_;
  float Mjj_;
  float ptJet1_;
  float ptJet2_;
  int cat_;

  float mPairPFCiC_;
  float mPairNoCorrPFCiC_;
  float mPairResPFCiC_;
  float mPairResWrongVtxPFCiC_;
  int diPhoVtxPFCiC_;
  float vtxProbPFCiC_;
  float diPhoVtxXPFCiC_;
  float diPhoVtxYPFCiC_;
  float diPhoVtxZPFCiC_;
  float MjjPFCiC_;
  float ptJet1PFCiC_;
  float ptJet2PFCiC_;
  int catPFCiC_;

  float mPairCiC_;
  float mPairNoCorrCiC_;
  float mPairResCiC_;
  float mPairResWrongVtxCiC_;
  int diPhoVtxCiC_;
  float vtxProbCiC_;
  float diPhoVtxXCiC_;
  float diPhoVtxYCiC_;
  float diPhoVtxZCiC_;
  float MjjCiC_;
  float ptJet1CiC_;
  float ptJet2CiC_;

  float cosThetaLead;
  float cosThetaLeadPFCiC;
  float cosThetaLeadCiC;
  
  Int_t nOutPhotons_;
  std::vector<ReducedPhotonData> OutPhotons_;
  Int_t nOutPhotonsPFCiC_;
  std::vector<ReducedPhotonData> OutPhotonsPFCiC_;
  Int_t nOutPhotonsCiC_;
  std::vector<ReducedPhotonData> OutPhotonsCiC_;

  float MET;
  float METPhi;

  VecbosPho pho1_;
  VecbosPho pho2_;

  int nSigma;
  std::vector<float> mPairScale;
  std::vector<float> pho1MVAScale;
  std::vector<float> pho2MVAScale;
  std::vector<float> diPhoMVAScale;
  std::vector<float> mPairSmear;
  std::vector<float> pho1MVASmear;
  std::vector<float> pho2MVASmear;
  std::vector<float> diPhoMVASmear;

  std::vector<float> mPairScalePFCiC;
  std::vector<float> mPairSmearPFCiC;

  std::vector<float> mPairScaleCiC;
  std::vector<float> mPairSmearCiC;

  //for mumuG
  const static int maxMuMuG = 500;
  int nMuMuG;
  float massMuMuGamma[maxMuMuG];
  float massMuMuRegGamma[maxMuMuG];
  float massMuMuScaleGamma[maxMuMuG];
  float massMuMuGenGamma[maxMuMuG];
  float massMuMu[maxMuMuG];
  float puWeight[maxMuMuG];
  
  MuCollection MMG_Mu1;
  MuCollection MMG_Mu2;
  PhoCollection MMG_Pho;
  float mvaPho[maxMuMuG];
  float isosumoetPho[maxMuMuG];
  
  //----
  float genHiggsPt;
  float genHiggsVx;
  float genHiggsVy;
  float genHiggsVz;

  float ptGenPho1;
  float etaGenPho1;
  float phiGenPho1;
  float energyGenPho1;
  float ptGenPho2;
  float etaGenPho2;
  float phiGenPho2;
  float energyGenPho2;
  
  float nPU_;
  int nVtxOut;
  int lumiBlockOut;
  int runNumberOut;
  int evtNumberOut;


  //member data
  const static int maxPho = 100;
  int nPho_;
  std::vector<VecbosPho> *Photons_; // this contains ALL photons
  
  bool photonMatchedElectron[maxPho];

  //this is the collection of the TMVA selected vertices for each photon pair
  //the format is std::pair< (iPho1 << 14) + iPho2,iVrt> 
  int nPair_;
  std::vector<std::pair<int,int> > *ggVerticesPhotonIndices;
  std::vector<std::pair<int, float> > *ggVerticesVertexIndex;
  std::vector<float>                  *ggVerticesPerEvtMVA;

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
