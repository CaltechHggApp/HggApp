#ifndef SusyHggSelector_h
#define SusyHggSelector_h

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

class SusyHggSelector : public BaseSelector {
public:
  SusyHggSelector(std::vector<std::string> fNames, std::string treeName,std::string outputFile);
  void setOptimize(){optimize=true;}
  void setIsMC(){isMC=true;}

  void setSmearingCFG(std::string smearCFG);

  void setIsDY(bool b = true){ isDY=b;}
protected:

  StandardPhotonID photonID;
  StandardElectronID electronID;
  VecbosJetID      jetID;

  JECUReader       jecReader;

  HggEnergyScale  *smearer=0;
  //mandatory overrides
  virtual void processEntry(Long64_t iEntry);
  virtual int init();
  virtual void processConfig(ReadConfig& cfg);
  virtual void setupOutputTree();

  virtual void clear();

  void selectJets(std::vector<TLorentzVector> *selectedJets,int correction); //+1 -> 1 sigma up JEC, -1 -> 1 sigma down JEC, 0 -> no JEC

  void fillGenTruth(); //fill the generator truth info in the output

  /*
  //functions to compute the BTagSF for the event
  int getptbin_for_btag(float pt);
  int get_eta_bin_jet(float eta);
  void get_SF_btag(const VecbosJet& jet, float &SF, float &SFerr);
  void computeBTagSF();
  */

  bool optimize=false; //setting this flag will reduce selection cuts and output more info
  bool isMC=false;
  bool isDY=false; // is a Z->ll sample

  //selection cuts
  float min_pho1_pt=32.;
  float min_pho2_pt=24.;
  float max_pho_eta=2.5;

  float min_jet_pt=30.;
  float max_jet_eta=3.;

  float min_mgg = 100;
  float max_mgg = 180;
  float min_ptgg = 20;

  //output variables
  int runNum;
  int lumiSec;
  int evtNum;

  float mgg;
  float ptgg;
  float etagg;
  float phigg;
  int hemgg;

  float pho1_pt;
  float pho1_eta;
  float pho1_phi;
  float pho1_r9;
  float pho1_seoe;
  bool pho1_genMatch;


  bool pho1_pass_id;
  bool pho1_pass_iso;
  bool pho1_pass_pixel;

  //optimize only
  float pho1_sieie;
  float pho1_HE;
  float pho1_charged;
  float pho1_neutral;
  float pho1_photon;
  bool pho1_eleveto;

  float pho1_energyGen;
  
  float pho2_pt;
  float pho2_eta;
  float pho2_phi;
  float pho2_r9;
  float pho2_seoe;
  bool pho2_genMatch;

  bool pho2_pass_id;
  bool pho2_pass_iso;
  bool pho2_pass_pixel;

  //optimize only
  float pho2_sieie;
  float pho2_HE;
  float pho2_charged;
  float pho2_neutral;
  float pho2_photon;
  bool pho2_eleveto;

  float pho2_energyGen;

  int nJ;
  int nJ_up;
  int nJ_down;
  
  float hem1_pt;
  float hem1_eta;
  float hem1_phi;
  float hem1_M;

  float hem2_pt;
  float hem2_eta;
  float hem2_phi;
  float hem2_M;

  //down 1 sigma JEC
  float hem1_pt_down;
  float hem1_eta_down;
  float hem1_phi_down;
  float hem1_M_down;

  float hem2_pt_down;
  float hem2_eta_down;
  float hem2_phi_down;
  float hem2_M_down;

  //up 1 sigma JEC
  float hem1_pt_up;
  float hem1_eta_up;
  float hem1_phi_up;
  float hem1_M_up;

  float hem2_pt_up;
  float hem2_eta_up;
  float hem2_phi_up;
  float hem2_M_up;


  float MET;
  float METphi;

  float HT;
  float HT_up;
  float HT_down;

  float MHT;
  float MHT_up;
  float MHT_down;

  float MHTphi;
  float MHTphi_up;
  float MHTphi_down;

  float t1MET;
  float t1METphi;

  float MR;
  float Rsq;
  float t1Rsq;
  
  float MR_down;
  float Rsq_down;
  float t1Rsq_down;
  
  float MR_up;
  float Rsq_up;
  float t1Rsq_up;
  
  float mu1_pt;
  float mu1_eta;
  float mu1_phi;

  float ele1_pt;
  float ele1_eta;
  float ele1_phi;

  float highest_csv;
  float highest_csv_pt;
  float highest_csv_eta;
  float highest_csv_phi;

  float highest_csv_down;
  float highest_csv_pt_down;
  float highest_csv_eta_down;
  float highest_csv_phi_down;

  float highest_csv_up;
  float highest_csv_pt_up;
  float highest_csv_eta_up;
  float highest_csv_phi_up;

  float second_csv;
  float second_csv_pt;
  float second_csv_eta;
  float second_csv_phi;

  float second_csv_down;
  float second_csv_pt_down;
  float second_csv_eta_down;
  float second_csv_phi_down;

  float second_csv_up;
  float second_csv_pt_up;
  float second_csv_eta_up;
  float second_csv_phi_up;

  float mbb;
  float mbb_up;
  float mbb_down;

  float mbb_NearH;
  float mbb_NearH_up;
  float mbb_NearH_down;

  float mbb_NearZ;
  float mbb_NearZ_up;
  float mbb_NearZ_down;

  float btagSF;
  float puWeight;

  const static size_t kMaxJets=20;
  int indexJet[kMaxJets];
  float ptJet[kMaxJets];
  float etaJet[kMaxJets];
  float phiJet[kMaxJets];
  float energyJet[kMaxJets];
  float corrUpJet[kMaxJets];
  float corrDownJet[kMaxJets];
  int hemJet[kMaxJets];

  const static int maxSusyPart=50;
  int nSusyPart;
  int idSusyPart[maxSusyPart];
  float mSusyPart[maxSusyPart];

  float m22=0;
  float m23=0;
  float m24=0;
  float m25=0;
};

#endif
