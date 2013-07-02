#include "include/McTruthEvent.hh"
#include <iostream>

using namespace std;

McTruthEvent::McTruthEvent() {
  cached_ = false;
  cachedTop_ = false;
  isData_ = false;
  idxEleZPrompt_ = -1;
  idxPosZPrompt_ = -1;
  idxEleWPrompt_ = -1;
  idxHadronTopPrompt_  = -1;
  idxEleTopPrompt_     = -1;
  idxPosTopPrompt_     = -1;
  idxMuTopPrompt_      = -1;
  idxAMuTopPrompt_     = -1;
  idxTauTopPrompt_     = -1;
  idxTauEleTopPrompt_  = -1;
  idxATauPosTopPrompt_ = -1;
  idxTauMuTopPrompt_   = -1;
  idxATauAMuTopPrompt_ = -1;
}

void McTruthEvent::LoadDecay(Int_t n, Int_t *id, Int_t *mother) {
  n_ = n;
  id_ = id;
  moth_ = mother;

  cached_ = false;
  cachedTop_ = false;
  idxEleZPrompt_ = -1;
  idxPosZPrompt_ = -1;
  idxEleWPrompt_ = -1;
  idxHadronTopPrompt_  = -1;
  idxEleTopPrompt_     = -1;
  idxPosTopPrompt_     = -1;
  idxMuTopPrompt_      = -1;
  idxAMuTopPrompt_     = -1;
  idxTauTopPrompt_     = -1;
  idxTauEleTopPrompt_  = -1;
  idxATauPosTopPrompt_ = -1;
  idxTauMuTopPrompt_   = -1;
  idxATauAMuTopPrompt_ = -1;
}

void McTruthEvent::LoadMomentum(Float_t *p, Float_t *energy, Float_t *theta, Float_t *phi) {
  p_ = p;
  energy_ = energy;
  theta_ = theta;
  phi_ = phi;
}

void McTruthEvent::SearchVecBos(int leptonPdgId) {

  // WARNING: loop over only the first 50 entries ot to use too much CPU
  for(int imc=0; imc<50; imc++) {
    if ( (id_[imc]==leptonPdgId)  && (fabs(id_[moth_[imc]])==23) ) idxEleZPrompt_=imc;
    if ( (id_[imc]==-1*leptonPdgId)  && (fabs(id_[moth_[imc]])==23) ) idxPosZPrompt_=imc;
    if ( (fabs(id_[imc])==leptonPdgId)  && (fabs(id_[moth_[imc]])==24) ) idxEleWPrompt_=imc;
  }
  cached_ = true;
}

void McTruthEvent::SearchTop() {

  for(int imc=0; imc<50; imc++) {
    int theid       = id_[imc];
    int themoth     = moth_[imc];
    int theGmoth    = moth_[themoth];
    int theGGmoth   = moth_[theGmoth];
    int themothid   = id_[themoth];
    int theGmothid  = id_[theGmoth];
    int theGGmothid = id_[theGGmoth];
    
    if ( (theid==11)  && 
	 ( ((abs(themothid)==24) && (abs(theGmothid)==6)) 
	   || 
	   ((abs(themothid)==24) && (abs(theGmothid)==24) && (abs(theGGmoth)==6)) 
	   ) ) idxEleTopPrompt_=imc;
    
    if ( (theid==-11)  && 
	 ( ((abs(themothid)==24) && (abs(theGmothid)==6)) 
	   || 
	   ((abs(themothid)==24) && (abs(theGmothid)==24) && (abs(theGGmoth)==6)) 
	   ) ) idxPosTopPrompt_=imc;

    if ( (theid==13)  && 
	 ( ((abs(themothid)==24) && (abs(theGmothid)==6)) 
	   || 
	   ((abs(themothid)==24) && (abs(theGmothid)==24) && (abs(theGGmoth)==6)) 
	   ) ) idxMuTopPrompt_=imc;

    if ( (theid==-13)  && 
	 ( ((abs(themothid)==24) && (abs(theGmothid)==6)) 
	   || 
	   ((abs(themothid)==24) && (abs(theGmothid)==24) && (abs(theGGmoth)==6)) 
	   ) ) idxAMuTopPrompt_=imc;

    if ( (theid==15 || theid==-15) && 
	 ( ((abs(themothid)==24) && (abs(theGmothid)==6))
	   || 
	   ((abs(themothid)==24) && (abs(theGmothid)==24) && (abs(theGGmoth)==6))
	   ) ) idxTauTopPrompt_=imc;

    if ( (abs(theid)==1 || abs(theid)==2 || abs(theid)==3 || abs(theid)==4)  && 
	 ( ((abs(themothid)==24) && (abs(theGmothid)==6)) 
	   || 
	   ((abs(themothid)==24) && (abs(theGmothid)==24) && (abs(theGGmoth)==6)) 
	   ) ) idxHadronTopPrompt_=imc;
  }

  // for tau decays, check if tau -> ele / mu or hadrons
  if (idxTauTopPrompt_>-1) {
    for(int imc=0; imc<700; imc++) {
      int theid     = id_[imc];
      int themoth   = moth_[imc];
      int themothid = id_[themoth];
      
      if ( theid== 11 && themothid== 15 ) idxTauEleTopPrompt_ =imc;
      if ( theid==-11 && themothid==-15 ) idxATauPosTopPrompt_=imc;
      if ( theid== 13 && themothid== 15 ) idxTauMuTopPrompt_  =imc;
      if ( theid==-13 && themothid==-15 ) idxATauAMuTopPrompt_=imc;
    }
  }

  cachedTop_ = true;
}

TLorentzVector McTruthEvent::p4ElectronWenuPrompt() {

  if(!cached_) SearchVecBos();

  if(idxEleWPrompt_ < 0) {
    TLorentzVector zero(0.,0.,0.,0.);
    return zero;
  }

  float energy = energy_[idxEleWPrompt_];
  TVector3 p3;
  p3.SetMagThetaPhi(p_[idxEleWPrompt_],theta_[idxEleWPrompt_],phi_[idxEleWPrompt_]);
  TLorentzVector v(p3,energy);

  return v;
}

std::pair<TLorentzVector,TLorentzVector> McTruthEvent::p4ElectronZeePrompt() {

  if(!cached_) SearchVecBos();

  if(idxEleZPrompt_ < 0 || idxPosZPrompt_ < 0) {
    TLorentzVector zero1(0.,0.,0.,0.);
    TLorentzVector zero2(0.,0.,0.,0.);
    return std::make_pair<TLorentzVector,TLorentzVector>(zero1,zero2);
  }

  float energy1 = energy_[idxEleZPrompt_];
  TVector3 p1;
  p1.SetMagThetaPhi(p_[idxEleZPrompt_],theta_[idxEleZPrompt_],phi_[idxEleZPrompt_]);
  TLorentzVector v1(p1,energy1);

  float energy2 = energy_[idxPosZPrompt_];
  TVector3 p2;
  p2.SetMagThetaPhi(p_[idxPosZPrompt_],theta_[idxPosZPrompt_],phi_[idxPosZPrompt_]);
  TLorentzVector v2(p2,energy2);
  
  return std::make_pair<TLorentzVector,TLorentzVector>(v1,v2);
}

int McTruthEvent::indexEleWPrompt(int leptonPdgId) {
  if(!cached_) SearchVecBos(leptonPdgId);
  return idxEleWPrompt_;
}

std::pair<int,int> McTruthEvent::indicesEleZPrompt(int leptonPdgId) {
  if(!cached_) SearchVecBos(leptonPdgId);
  return std::make_pair<int,int>(idxEleZPrompt_,idxPosZPrompt_);
}

int McTruthEvent::indexEleTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxEleTopPrompt_;
}

int McTruthEvent::indexPosTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxPosTopPrompt_;
}

int McTruthEvent::indexMuTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxMuTopPrompt_;
}

int McTruthEvent::indexAMuTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxAMuTopPrompt_;
}

int McTruthEvent::indexTauEleTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxTauEleTopPrompt_;
}

int McTruthEvent::indexATauPosTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxATauPosTopPrompt_;
}

int McTruthEvent::indexTauMuTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxTauMuTopPrompt_;
}

int McTruthEvent::indexATauAMuTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxATauAMuTopPrompt_;
}

int McTruthEvent::indexHadronTopPrompt() {
  if(!cachedTop_) SearchTop();
  return idxHadronTopPrompt_;
}
