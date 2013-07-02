//--------------------------------------------------------------
// Description:
//    Class for contain the ee MC truth informations
// Authors:
//    Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//--------------------------------------------------------------

#ifndef MCTRUTHEVENT_H
#define MCTRUTHEVENT_H

#include <TROOT.h>
#include <TLorentzVector.h>


class McTruthEvent {

public:
  McTruthEvent();
  ~McTruthEvent() {};

  void SetData(bool what) { isData_ = what; };
  void LoadDecay(Int_t n, Int_t *id, Int_t *mother);
  void LoadMomentum(Float_t *p, Float_t *energy, Float_t *theta, Float_t *phi);

  TLorentzVector p4ElectronWenuPrompt();
  std::pair<TLorentzVector,TLorentzVector> p4ElectronZeePrompt();

  std::pair<int,int> indicesEleZPrompt(int leptonPdgId=11);
  int indexEleWPrompt(int leptonPdgId=11);
  int indexEleTopPrompt();
  int indexPosTopPrompt();
  int indexMuTopPrompt();
  int indexAMuTopPrompt();
  int indexTauTopPrompt();
  int indexATauTopPrompt();
  int indexTauEleTopPrompt();
  int indexATauPosTopPrompt();
  int indexTauMuTopPrompt();
  int indexATauAMuTopPrompt();
  int indexHadronTopPrompt();

private:

  void SearchVecBos(int leptonPdgId=11);
  void SearchTop();

  bool cached_;
  bool cachedTop_;
  bool isData_;

  Int_t n_;
  Int_t *id_;
  Int_t *moth_;
  
  Float_t *p_, *energy_, *theta_, *phi_;
  int idxEleZPrompt_, idxPosZPrompt_, idxEleWPrompt_;
  int idxEleTopPrompt_, idxPosTopPrompt_;
  int idxMuTopPrompt_, idxAMuTopPrompt_;
  int idxHadronTopPrompt_;
  int idxTauTopPrompt_;
  int idxTauEleTopPrompt_, idxATauPosTopPrompt_; 
  int idxTauMuTopPrompt_, idxATauAMuTopPrompt_; 

};

#endif

