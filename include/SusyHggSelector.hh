#ifndef SusyHggSelector_h
#define SusyHggSelector_h

#include "BaseSelector.hh"
#include "StandardPhotonID.hh"
#include "StandardElectronID.hh"
#include "VecbosJetID.hh"
#include "RazorVariables.hh"

#include "TVector3.h"
#include "TLorentzVector.h"

#include "assert.h"
#include <bitset>

class SusyHggSelector : public BaseSelector {
public:
  SusyHggSelector(std::vector<std::string> fNames, std::string treeName,std::string outputFile):BaseSelector(fNames,treeName,outputFile){setDoFill(false);}
  void setOptimize(){optimize=true;}
  void setIsMC(){isMC=true;}
protected:

  StandardPhotonID photonID;
  StandardElectronID electronID;
  VecbosJetID      jetID;
  //mandatory overrides
  virtual void processEntry(Long64_t iEntry);
  virtual int init();
  virtual void processConfig(ReadConfig& cfg);
  virtual void setupOutputTree();

  virtual void clear();

  void fillGenTruth(); //fill the generator truth info in the output

  bool optimize=false; //setting this flag will reduce selection cuts and output more info
  bool isMC=false;

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
  float mgg;
  float ptgg;
  float etagg;
  float phigg;

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

  int nJ;
  
  float hem1_pt;
  float hem1_eta;
  float hem1_phi;
  float hem1_M;

  float hem2_pt;
  float hem2_eta;
  float hem2_phi;
  float hem2_M;

  float MET;
  float METphi;

  float MR;
  float Rsq;
  
  float mu1_pt;
  float mu1_eta;
  float mu1_phi;

  float ele1_pt;
  float ele1_eta;
  float ele1_phi;

  float puWeight;

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
