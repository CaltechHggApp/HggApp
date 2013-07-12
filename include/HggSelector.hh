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

#include "BaseSelector.hh"

struct OutputVars{ // all the misc vars to write to the output tree for each selection
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
  float cosThetaLead_;

  int nOutPhotons_;
  std::vector<ReducedPhotonData> OutPhotons_;

  std::vector<float> mPairScale_;
  std::vector<float> mPairSmear_;

  void clear(){ 
    OutPhotons_.clear(); mPairScale_.clear(); mPairScale_.clear();
    mPair_=0; mPairNoCorr_=0; mPairRes_=0; mPairResWrongVtx_=0;
    diPhoMVA_=-999; diPhoVtx_=-1; vtxProb_=0;
    diPhoVtxX_=-999; diPhoVtxY_=-999; diPhoVtxZ_=-999;
    Mjj_=-1; ptJet1_=-1; ptJet2_=-1;
    cat_=-1; cosThetaLead_=-2;
  }
}; 

class HggSelector: public BaseSelector{
public:
  HggSelector():BaseSelector(){}
  HggSelector(vector<string> fNames,string treeName,string outputFile):
    BaseSelector(fNames,treeName,outputFile){}

  static constexpr unsigned int nSelections=3;
  enum SelectionType: unsigned int{kMVA=0,kPFCiC=1,kCiC=2};

  void suppressElectronVeto(){doElectronVeto=false;}
  void setForceVertexZero(){forceVtxZero=true;}

protected:
  //methods we need/want to override
  virtual int init() final;
  virtual void clear() final;
  virtual void setDefaults() final;
  virtual void setupOutputTree() final;
  virtual void processConfig(ReadConfig& cfg) final;
  virtual void firstInit() final;

  virtual void processEntry(Long64_t iEntry) final; //main override
  void processOnce(); // only run this part once, regardless of how many selections are turned on
  void processEntry(SelectionType SelType); // run this once per selection


  //Hgg Specific methods
  void setMassResConfig(string s){massResConfig = s;}
  void setNSigma(int ns){nSigma=ns;}

  bool preSelectPhotons(VecbosPho*,VecbosPho*,TVector3); // kinematic photons selections
  void smearPhoton(VecbosPho*,int);

  std::pair<int,int> getBestPair(float*,int,int);
  std::pair<int,int> getBestPairCiC(int,int,bool);
  int getCategory(OutputVars & vars);
  int getCategoryPFCiC(OutputVars & vars);

  void setupTMVA();
  void fillGenInfo();

  float getVBFMjj(VecbosPho*,VecbosPho*,TVector3,float*);
  bool passJetID(VecbosJet*);

  //get vertex info
  int getVertexIndex(int,int);
  float getVertexProb(int,int);
  float getVertexMVA(int,int);

  float getDiPhoMVA(int,int,float,float,bool);

  float getMPair(int,int);

  ReducedPhotonData getReducedData(VecbosPho*,TVector3,int);


  //Hgg specific flags
  bool doElectronVeto = true;
  bool forceVtxZero   = false;

  //we can re-run the regression, scaling or smearing
  bool doRegression;
  bool doScale;
  bool doSmear;
  HggEGEnergyCorrector* corrector;
  HggEnergyScale* scale;
  int applyScaleSmear;
  HggEnergyScale* smear;
  
  //Mass Resolution
  string massResConfig;
  HggMassResolution *massRes;

  //Photon ID
  HggPhotonID *PhotonID;

  //selection
  float leadPhoEtMin;
  float subleadPhoEtMin;
  bool doPtOverM;
  float PtOverMLead;
  float PtOverMSubLead;
  static constexpr float rhoFac    = 0.17;
  static constexpr float rhoFacBad = 0.52;
  static constexpr float isoSumConst = 0;//5;
  static constexpr float isoSumConstBad = 0;//7;

  //TMVA stuff
  string weightFile_diPho;
  string methodName_diPho;
  TMVA::Reader *diPhotonMVA;

  std::map<std::string,TH1F*> MVAInputs;
  std::map<std::string,TH1F*> MjjDists;
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

  //outputSets:
  OutputVars vars[nSelections];

  //Output Variables:
  int trigger_;

  /*
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
*/
  float MET;
  float METPhi;

  VecbosPho pho1_;
  VecbosPho pho2_;

  int nSigma = 3; 

  std::vector<float> pho1MVAScale;
  std::vector<float> pho2MVAScale;
  std::vector<float> diPhoMVAScale;

  std::vector<float> pho1MVASmear;
  std::vector<float> pho2MVASmear;
  std::vector<float> diPhoMVASmear;

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
};
