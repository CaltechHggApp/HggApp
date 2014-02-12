#include "VecbosJetID.hh"
#include <iostream>

VecbosJetID::VecbosJetID() {
  eFracCut[kLoose]  = 0.99;
  eFracCut[kMedium] = 0.95;
  eFracCut[kTight]  = 0.90;
}

bool VecbosJetID::passID(const VecbosJet& jet, WP wp) {
  /*
  if( jet.chargedHadronMultiplicity==0 &&
      jet.neutralHadronMultiplicity==0 &&
      jet.photonMultiplicity==0 &&
      jet.electronMultiplicity==0 &&
      jet.muonMultiplicity==0 &&
      jet.HFHadronMultiplicity==0 &&
      jet.HFEMMultiplicity==0) return false;

  */

  //if(jet.photonEnergy/jet.energy > eFracCut[wp]) return false; //neutral EM energy
  if(jet.neutralHadronFraction > eFracCut[wp]) return false;
  if(fabs(jet.eta) < 2.4) {
    //if(jet.chargedHadronMultiplicity == 0) return false;
    if(jet.chargedHadronFraction == 0) return false;
    //if(jet.electronEnergy/jet.energy > 0.99) return false;
  }

  return true;
}
