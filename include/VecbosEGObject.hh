#ifndef VecbosEGObject_h
#define VecbosEGObject_h

#include "VecbosBase.hh"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "../src/HggPhysUtils.cc"

//object for a basic cluster
class VecbosBC{
public:
  VecbosBC();
  VecbosBC(VecbosBase*, int);
  virtual void Init(VecbosBase*,int);
  int index;
  float energy;
  float eta;
  float phi;

  float etaCrystal;
  float phiCrystal;
  int iEta;
  int iPhi;
  float thetaTilt;
  float phiTilt;

};

class VecbosPFBC : public VecbosBC{
public:
  VecbosPFBC();
  VecbosPFBC(VecbosBase*, int);
  void Init(VecbosBase* o, int i);
};

//super cluster
class VecbosSC{
public:
  VecbosSC();
  VecbosSC(VecbosBase*, int);
  virtual void Init(VecbosBase*, int);
  VecbosBC BCSeed;
  int nBCs;
  std::vector<VecbosBC> basicClusters; // basic clusters associated with this supercluser
  int index;
  float energy;
  float esEnergy;
  float eta;
  float phi;
  float e3x3;
  float e5x5;

  float e3x1;
  float e1x3;
  float e4x4;
  float eMax;
  float e2x2;
  float e2nd;
  float e1x5;
  float e2x5Max;
  float e2x5Left;
  float e2x5Right;
  float e2x5Top;
  float e2x5Bottom;
  
  float eLeft;
  float eRight;
  float eTop;
  float eBottom;

  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;

  float esEffSigRR;
  
  TVector3 CaloPos;

  float rawE;
  float phiWidth;
  float etaWidth;
  float HoverE;

  //for Hgg corrector

  float r9;
  float s4ratio;
  //  SCInfo getStruct();
};

class VecbosPFSC : public VecbosSC{
public:
  VecbosPFSC();
  VecbosPFSC(VecbosBase*, float, float); // for eta/phi match to the SC position
  VecbosPFSC(VecbosBase*, int);
  void Init(VecbosBase*, int);
  void Init(VecbosBase*, float, float);
  std::vector<VecbosPFBC> pfClusters;
};

class VecbosConversion{
public:
  VecbosConversion();
  VecbosConversion(VecbosBase*,int);
  void Init(VecbosBase*, int);
  int index;
  TVector3 pPair;
  TVector3 pRefittedPair;
  TLorentzVector p4RefittedPair;

  TVector3 CaloPos;

  float eOverP;

  TVector3 vtx;
  float vtxChi2;
  float vtxChi2Prob;
  bool vtxIsValid;
  int vtxNTracks;
  float vtxMVA;

  //track info
  float trk1Dz;       
  float trk1DzError;  
  float trk1Charge;   
  float trk1Algo;     
  float trk1D0;       
  float trk1Pout;
  float trk1Pin;

  float trk2Dz;       
  float trk2DzError;  
  float trk2Charge;   
  float trk2Algo;     
  float trk2D0;       
  float trk2Pout;     
  float trk2Pin; 

  //ConversionInfo getStruct();
};

class VecbosGen{
public:
  VecbosGen();
  VecbosGen(VecbosBase*, int);
  void Init(VecbosBase*, int);
  int index;
  float energy;
  float pt;
  float eta;
  float phi;
  float mass;
  //TLorentzVector getP4();

  //TVector3 Vtx;
  float Vx;
  float Vy;
  float Vz;

  int status;
  int id;
  int statusMother;
  int idMother;
  int indexMother;
};


class VecbosPho{
public:
  VecbosPho();
  VecbosPho(VecbosBase*,int);
  void Init(VecbosBase*, int);
  void matchConversion(VecbosBase*,bool);
  int index;
  float energy;
  float eta;
  float phi;


  //Hgg Correction Variables
  float correctedEnergy;
  float correctedEnergyError;
  float scaledEnergy;
  float scaledEnergyError;
  //float smearedEnergy;
  //float smearedEnergyError;
  float finalEnergy;
  float finalEnergyError;

  float dEoE;
  float dEoEErr;
  VecbosSC SC;
  VecbosPFSC PFSC;
  float HoverE;
  float HTowOverE;

  int hasPixel;
  TVector3 CaloPos;

  TLorentzVector p4FromVtx(TVector3 vtx,float E,bool pf=false);
  VecbosConversion conversion;
  //isolation variables
  float dr03EcalRecHitSumEtCone;
  float dr03HcalTowerSumEtCone;
  float dr03TrkSumPtCone;
  float dr03TrkSumPtHollowCone;
  float dr04EcalRecHitSumEtCone;
  float dr04HcalTowerSumEtCone;
  float dr04TrkSumPtCone;
  float dr04TrkSumPtHollowCone;

  int           nPV;
  //tracker isolation
  float        dr02TrackIso[100];
  float        dr03TrackIso[100];
  float        dr04TrackIso[100];

  //pfIsolation
  float         dr01ChargedHadronPFIso[100];
  float         dr02ChargedHadronPFIso[100];
  float         dr03ChargedHadronPFIso[100];
  float         dr04ChargedHadronPFIso[100];
  float         dr05ChargedHadronPFIso[100];
  float         dr06ChargedHadronPFIso[100];

  float                      dr01NeutralHadronPFIso;
  float                      dr02NeutralHadronPFIso;
  float                      dr03NeutralHadronPFIso;
  float                      dr04NeutralHadronPFIso;
  float                      dr05NeutralHadronPFIso;
  float                      dr06NeutralHadronPFIso;
  float                      dr01PhotonPFIso;
  float                      dr02PhotonPFIso;
  float                      dr03PhotonPFIso;
  float                      dr04PhotonPFIso;
  float                      dr05PhotonPFIso;
  float                      dr06PhotonPFIso;

  bool isBarrel(){return (fabs(this->SC.eta) < 1.48);}
  int  getCategory(){ return (SC.r9>0.94)+2*(isBarrel()); } //get the category 0-3 of the photon
  //PhoInfo getStruct();

  VecbosGen genMatch;
  void doGenMatch(VecbosBase* o);
};

struct ReducedPhotonData{
  Float_t pt,eta,phi,E,EError,EErrorSmeared;
  Float_t pt_Gen, eta_Gen, phi_Gen, E_Gen;
  Float_t etaSC;
  int index;
  float r9;
  bool passPFCiC;
  int category;
  float idMVA;
  int mother;
  float HoverE;
  float sieie;
  float dr03PFChargedIso;
  float isosumGood;
  float isosumBad;
  float dr03EcalIso;
  float dr04HcalIso;
  float dr03TrackIso;
  float dr02PFChargedIso;

};

class VecbosEle{
public:
  VecbosEle();
  VecbosEle(VecbosBase*,int);
  void Init(VecbosBase*, int);
  int index;
  float pt;
  float energy;
  float eta;
  float phi;

  TLorentzVector getP4(float E){
    TLorentzVector p4;
    p4.SetPtEtaPhiM(E/cosh(eta),eta,phi,0);
    return p4;
  }
  int charge;

  float correctedEnergy;
  float correctedEnergyError;

  VecbosSC SC;
  float esEnergy;
  float HoverE;
  bool isEcalDriven;
  bool isTrackerDriven;
  
  float vtxX;
  float vtxY;
  float vtxZ;

  float EOverP;

  float d0Track;
  float dzTrack;

  float dEtaSCTrackAtVtx;
  float dPhiSCTrackAtVtx;

  float dEtaSCTrackAtCalo;
  float dPhiSCTrackAtCalo;

  float dr03ChargedHadronPFIso;
  float dr03NeutralHadronPFIso;
  float dr03PhotonPFIso;

  float dr04ChargedHadronPFIso;
  float dr04NeutralHadronPFIso;
  float dr04PhotonPFIso;

  float dr03TkSumPt;
  float dr03EcalRecHitSumEt;
  float dr03HcalTowerSumEt;

  float dr04TkSumPt;
  float dr04EcalRecHitSumEt;
  float dr04HcalTowerSumEt;

  float idMVA;

  bool hasMatchedConversion;
  int expInnerLayersHits;

  VecbosGen genMatch;
  void doGenMatch(VecbosBase*);

  //EleInfo getStruct();
};

class VecbosMu{
public:
  VecbosMu();
  VecbosMu(VecbosBase*, int);
  void Init(VecbosBase*, int);
  int index;
  float energy;
  float pt;
  float eta;
  float phi;
  //TLorentzVector p4;
  int charge;
  float combinedIso;
  float emIso;
  float hadIso;
  float trkIso;
  bool isGlobalMuon;
  bool isTrackerMuon;
  bool isPromptMuon;
  int nTrackHits;
  int nPixelHits;
  float trackImpactPar;

  bool isLooseMuon;
  bool isTightMuon;

  VecbosGen genMatch;
  void doGenMatch(VecbosBase*);
};

class VecbosJet{
public:
  enum JetType{PFPUcorr,PFNoPU};
  VecbosJet();
  VecbosJet(VecbosBase*, int,JetType);
  void Init(VecbosBase*, int,JetType);
  int index;
  float energy;
  float uncorrEnergy;
  float uncorrpx;
  float uncorrpy;
  float uncorrpz;
  float pt;
  float eta;
  float phi;
  
  TLorentzVector getP4();

  int charge;

  VecbosJet::JetType type;
  
  float vtxX;
  float vtxY;
  float vtxZ;
  TVector3 getVertex(){return TVector3(vtxX,vtxY,vtxZ);}

  float area;
  float chargedHadronFraction;
  float neutralHadronFraction;

  float jetIdMva;
  
  float betaStar;
  float betaStarIdMVA;
  float betaStarClassicIdMVA;
  float rmsCands;
  float rmsCandsHand;

  float combinedSecondaryVertex;
  float simpleSecondaryVertexHighPur;
  float simpleSecondaryVertexHighEff;
};

typedef std::vector<VecbosPho> PhoCollection;
typedef std::vector<VecbosSC> SCCollection;
typedef std::vector<VecbosPFSC> PFSCCollection;
typedef std::vector<VecbosBC> BCCollection;
typedef std::vector<VecbosConversion> ConvCollection;
typedef std::vector<VecbosMu> MuCollection;
typedef std::vector<VecbosEle> EleCollection;
typedef std::vector<VecbosGen> GenCollection;
typedef std::vector<VecbosJet> JetCollection;
#endif
