//--------------------------------------------------------------
// Description:
//    Class for clean and count jets
// Authors:
//    Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//--------------------------------------------------------------

#ifndef JETCOUNTER_H
#define JETCOUNTER_H

#include <vector>
#include <TVector3.h>
#include "include/Jet.hh"
#include "include/BTagJet.hh"

class JetCounter {

public:
  JetCounter(std::vector<Jet> theJets);
  ~JetCounter() {};
  void SetBTagJets(std::vector<BTagJet> btagjets);
  void SetParticlesToRemove(std::vector<TVector3> particles) {particlesToBeRemoved_=particles;} ;
  void SetThresholds(float etThreshold, float absEtaThreshold);
  void SetDistance(float deltar) {drMin_=deltar;} ;
  void SetMinHadFrac(float frac) {hFracMin_=frac;} ;
  void SetPFJetId(int WP) {pfJetId_=WP;} ;
  void SetVerbose(bool what) {verbose_=what;} ;
  std::vector<Jet> getGoodJets();
  std::vector<BTagJet> getGoodBTagJets();
  int numGoodJets();
  int numBTaggedJets(int btagalgo, float val);
  std::vector<Jet> getGoodBTaggedJets(int btagalgo, float val);
  TVector3 sumPtGoodJets();
  void reset() { cached_ = false; }

private:
  float etThr_, absEtaThr_;
  float drMin_;
  float hFracMin_;
  int pfJetId_;

  //! input jets and b-tagging blocks
  std::vector<Jet> jets_;
  std::vector<BTagJet> btagjets_;

  //! output of cleaned jets and b-tagging blocks
  std::vector<Jet> goodJets_;
  std::vector<BTagJet> goodBTagJets_;

  //! output of b-tagged jets with a certain algorithm (see numBTaggedJets function)
  std::vector<Jet> goodBTaggedJets_;

  std::vector<TVector3> particlesToBeRemoved_;
  TVector3 sumPt_;

  bool cached_;
  bool verbose_;

};

#endif

