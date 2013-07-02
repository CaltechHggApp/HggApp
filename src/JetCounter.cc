#include <iostream>
#include "CommonTools/include/BTagAlgoBits.h"
#include "include/JetCounter.hh"

using namespace std;

JetCounter::JetCounter(std::vector<Jet> theJets) {

  drMin_ = 0.3;
  etThr_ = 0.0;
  absEtaThr_ = 1000.;
  hFracMin_ = 0.0;
  pfJetId_ = Jet::none;

  jets_ = theJets;
  sumPt_.SetXYZ(0.,0.,0.);
  
  cached_ = false;
  verbose_ = false;

}

void JetCounter::SetBTagJets(std::vector<BTagJet> btagjets) {
  btagjets_ = btagjets;
  if(btagjets_.size() != jets_.size()) cout << "WARNING: nasty error: the size of the b-tagging collection"
                                            << " is different from the jet collection size." << endl;
}

std::vector<Jet> JetCounter::getGoodJets() {

  if(cached_) return goodJets_;

  goodJets_.clear();
  goodBTagJets_.clear();
  sumPt_.SetXYZ(0.,0.,0.);

  if(verbose_) {
    cout << "===> Loop over jets starts..." << endl;
    cout << "===> " << particlesToBeRemoved_.size() << " particles to be removed:" << endl;
    std::vector<TVector3>::iterator particle;
    for(particle=particlesToBeRemoved_.begin(); particle!=particlesToBeRemoved_.end(); particle++) {
      particle->Print();
    }
  }

  int j=-1;
  std::vector<Jet>::iterator jet;
  for(jet=jets_.begin(); jet!=jets_.end(); jet++) {
    j++;
    TVector3 p3Jet = jet->Get3Vector(); 

    if(verbose_) {
      cout << "Jet candidate:" << endl;
      p3Jet.Print();
    }

    bool goodJet = true;

    // cleaning class: check if the particles to remove fall into the jet  
    std::vector<TVector3>::iterator particle;
    for(particle=particlesToBeRemoved_.begin(); particle!=particlesToBeRemoved_.end(); particle++) {
      float deltar = p3Jet.DeltaR(*particle);

      if(verbose_) cout << "\tthis jet has deltaR = " << deltar << endl;

      if(deltar < drMin_) {
        goodJet = false;
        break;
      }
    }
    
    if( !goodJet ) continue;

    if(verbose_) cout << "\tnot removed yet. Now checking |eta|<" << absEtaThr_ << "  and eT>" << etThr_ 
                      << " GeV  and hadfrac > " << hFracMin_ << endl;

    // acceptance
    if( fabs(p3Jet.Eta()) > absEtaThr_ ) continue;
    if( p3Jet.Pt() < etThr_ ) continue;

    // cut on hadronic fraction
    // if( jet->HadFrac() < hFracMin_ ) continue;

    // PF jet ID
    if(! jet->isPFJetID(pfJetId_) ) continue;

    if(verbose_) cout << "This is a good jet." << endl;
    
    // this is a good jet
    goodJets_.push_back(*jet); 
    if (btagjets_.size() > 0) goodBTagJets_.push_back(btagjets_[j]);
   sumPt_ += p3Jet.Pt();

  }  

  // to not recalculate it, use the cached values
  cached_ = true;
  return goodJets_;

}

std::vector<BTagJet> JetCounter::getGoodBTagJets() {
  getGoodJets();
  return goodBTagJets_;
}

void JetCounter::SetThresholds(float etThreshold, float absEtaThreshold) {
  etThr_=etThreshold;
  absEtaThr_=absEtaThreshold;
}

int JetCounter::numGoodJets() {
  std::vector<Jet> goodjets = getGoodJets();
  return goodjets.size();
}

TVector3 JetCounter::sumPtGoodJets() {
  std::vector<Jet> goodjets = getGoodJets();
  return sumPt_;
}

int JetCounter::numBTaggedJets(int btagalgo, float val) {
  int num=0;
  goodBTaggedJets_.clear();
  for(unsigned int j=0; j<goodBTagJets_.size(); ++j) {
    float btag;
    switch (btagalgo) {
    case bits::combinedSecondaryVertexBJetTags: 
      btag = goodBTagJets_[j].combinedSecondaryVertexBJetTags;
      break;
    case bits::combinedSecondaryVertexMVABJetTags:
      btag = goodBTagJets_[j].combinedSecondaryVertexMVABJetTags;
      break;
    case bits::jetBProbabilityBJetTags:
      btag = goodBTagJets_[j].jetBProbabilityBJetTags;
      break;
    case bits::simpleSecondaryVertexBJetTags:
      btag = goodBTagJets_[j].simpleSecondaryVertexBJetTags;
      break;
    case bits::softMuonBJetTags:
      btag = goodBTagJets_[j].softMuonBJetTags;
      break;
    case bits::trackCountingHighPurBJetTags:
      btag = goodBTagJets_[j].trackCountingHighPurBJetTags;
      break;
    case bits::trackCountingHighEffBJetTags:
      btag = goodBTagJets_[j].trackCountingHighEffBJetTags;
      break;
    default:
      cout << "ERROR: BTag algorithm: " << btagalgo << " is not recognized. Check code in: JetCounter.hh" << endl;
      btag = -10000.;
    }
    if(btag >= val) {
      goodBTaggedJets_.push_back(goodJets_[j]);
      num++;
    }
  }
  return num;
}

std::vector<Jet> JetCounter::getGoodBTaggedJets(int btagalgo, float val) {
  int numBJ = numBTaggedJets(btagalgo,val);
  return goodBTaggedJets_;
}
