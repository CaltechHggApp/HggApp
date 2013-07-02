// VecbosApp includes
#include "include/Jet.hh"
#include <iostream>

Jet::Jet(TLorentzVector p, double efrac, double hfrac){
  v = p;
  Efrac = efrac;
  Hfrac = hfrac;
  
  // default jet ID values
  neutralHadFrac = 0.0;
  neutralEmFraction = 0.0;
  nConstituents = 0;
  chargedHadFraction = 0.0;
  chargedMultiplicity = 0.0;
  chargedEmFraction = 0.0;  
}

Jet::~Jet() {}

void Jet::setPFJetID(double nHFrac, double nEmFrac, int nConst, double chHFrac, double chMult, double chEmFrac) {
  neutralHadFrac = nHFrac;
  neutralEmFraction = nEmFrac;
  nConstituents = nConst;
  chargedHadFraction = chHFrac;
  chargedMultiplicity = chMult;
  chargedEmFraction = chEmFrac;
}

bool Jet::isPFJetID(int WP) {
  switch(WP) {
  case none:
    return true;
    break;
  case loose:
    if(neutralHadFrac>=0.99 || neutralEmFraction>=0.99 || nConstituents<=1) return false;
    if(fabs(v.Eta())<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  case medium:
    if(neutralHadFrac>=0.95 || neutralEmFraction>=0.95 || nConstituents<=1) return false;
    if(fabs(v.Eta())<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  case tight:
    if(neutralHadFrac>=0.90 || neutralEmFraction>=0.90 || nConstituents<=1) return false;
    if(fabs(v.Eta())<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  default:
    std::cout << "Jet::isPFJetID(nt WP). Requested wrong Working point. Available are loose, medium, tight." << std::endl;
    return false;
  }
  return true;
}
