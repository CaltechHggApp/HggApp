 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef HggReducer_h
#define HggReducer_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

#include "TH1F.h"
#include "TFile.h"

#include "HggVertexingNew.hh"
#include "HggEGEnergyCorrector.hh"
#include "HggEnergyScale.hh"
#include "HggScaling.hh"
#include "VecbosJetCorrector.hh"

#include <exception>
#include <stdexcept>
using namespace std;

class HggReducer : public Vecbos{
public:
  HggReducer(TTree *tree=0); /// Class Constructor
  HggReducer(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false,int mod=-1); /// Class Constructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void addTrigger(string s){triggerNames.push_back(s);}
  void setConfig(string s){config = s;}
  void setCorrectionType(int s){correctionType = s;}
  void setScaleSmear(int s){applyScaleSmear = s;}
  void SetVtxConfigFile(string s){vertexCFG = s;}
  void SetScaleConfigFile(string s){energyScaleCFG = s;}
  void SetSmearConfigFile(string s){energySmearCFG = s;}
  void SetMinPhoSelection(int n){minPhoSel = n;}
private:
  bool _goodRunLS; 
  bool _isData;
  float _weight;
  TTree *_treeCond;

  float minPhoSel;
  //configuration file and MVA parameters
  string config;
  string vertexCFG;
  string energyScaleCFG;
  string energySmearCFG;
  string mcScalingCFG;

  void init(string); // do variable initialization 
  void clearAll();
  void setOutputBranches();
  
  float computeTrackIso(int iPho, int iVtx,
			float ptMinTrkIso,
			float outerConeTrkIso,
			float innerConeTrkIso,
			float etaStripHalfWidthTrkIso,
			float dzMaxTrkIso,
			float dyMaxTrkIso);

  //photon preselection variables
  void fillMuons();
  void fillElectrons();
  void fillJets();

  HggVertexing *vertexer;
  //energy correction variables
  string correctionType;
  HggEGEnergyCorrector *corrector;
  HggEGEnergyCorrector *elecorrector;

  //jet corrector
  VecbosJetCorrector *jetCorr;
  bool correctJets;

  //energy smearing
  int applyScaleSmear;
  HggEnergyScale *energyScale;
  HggEnergyScale *energySmear;

  HggScaling *scaler;
  TTree * outTree;
  // define variables for the output tree:
  vector<string> triggerNames; // list of all the triggers to consider
  int * triggerBits;       // this will be an array of the trigger decision per event (for the output tree)

  // ...
  //Event info
  int lumiBlockO; 
  int runNumberO; 
  int evtNumberO; 
  int bunchX; 
  int orbitNumber; 
  int evtTime; 

  int phyDeclared; 
  float rho; 
  float rhoEtaMax44; 

  vector<short> *pileupBunchX; 
  vector<short> *pileupNInteraction; 
  float pileupTrueNumInterations;

  //main object collections for the reduced tree 
  const static int maxPho=100;
  int nPho_;
  vector<VecbosPho>  Photons_;

  void matchPhotonsElectrons();
  bool photonMatchedElectron[maxPho];

  int nMu_;
  MuCollection Muons_;
  int nEle_;
  EleCollection Electrons_;
  int nJet_;
  JetCollection Jets_;

  //this is the collection of the TMVA selected vertices for each photon pair
  //the format is std::pair< (iPho1 << 14) + iPho2,iVrt> 
  int nPair_;
  vector<pair<int,int> > ggVerticesPhotonIndices; 
  vector<pair<int,float> > ggVerticesVertexIndex01;  // vertex, MVA score
  vector<pair<int,float> > ggVerticesVertexIndex02;  // vertex, MVA score
  vector<pair<int,float> > ggVerticesVertexIndex03;  // vertex, MVA score
  vector<float>            ggVerticesPerEvtMVA;
  //std::vector<int> ggVertices_

  //vertex information
  void fillVertexInfo();
  int nVtx; 
  static const int maxVtx = 100;
  float vtxX[maxVtx];
  float vtxY[maxVtx];
  float vtxZ[maxVtx];
  float vtxChi2[maxVtx];
  float vtxNdof[maxVtx];
  float vtxNormalizedChi2[maxVtx];
  int vtxTrackSize[maxVtx];
  int vtxIsFake[maxVtx];
  int vtxIsValid[maxVtx];


  //for pileup reweighting
  TH1F* pileupWeightHist;
  TFile *pileupWeightFile;
  float pileupWeight;

  //GENERATOR information
  virtual void fillGeneratorInfo();  
  static const int MAXGenSaved = 1000;
  //gen-leve phton
  int nGenPho;
  GenCollection GenPhotons;

  int nGenMu;
  GenCollection GenMuons;

  int nGenEle;
  GenCollection GenElectrons;

  int nGenHiggs;
  GenCollection GenHiggs;

  int nGenOthers;
  GenCollection GenOthers;

  int procID;
  float qScale;
  float nPu;

  float caloMet;
  float caloMetPhi;

  float pfMet;
  float pfMetPhi;

  float tcMet;
  float tcMetPhi;

  float pfMetType1;
  float pfMetType1Phi;

  float ECALLaserFilter;
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
#endif


