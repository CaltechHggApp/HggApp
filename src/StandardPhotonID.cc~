#include "StandardPhotonID.hh"

StandardPhotonID::StandardPhotonID() {

  cut_HE[std::make_pair(kLoose,kEB)] = 0.05;
  cut_HE[std::make_pair(kMedium,kEB)] = 0.05;
  cut_HE[std::make_pair(kTight,kEB)] = 0.05;

  cut_HE[std::make_pair(kLoose,kEE)] = 0.05;
  cut_HE[std::make_pair(kMedium,kEE)] = 0.05;
  cut_HE[std::make_pair(kTight,kEE)] = 0.05;

  cut_sieie[std::make_pair(kLoose,kEB)] = 0.012;
  cut_sieie[std::make_pair(kMedium,kEB)] = 0.011;
  cut_sieie[std::make_pair(kTight,kEB)] = 0.011;

  cut_sieie[std::make_pair(kLoose,kEE)] = 0.034;
  cut_sieie[std::make_pair(kMedium,kEE)] = 0.033;
  cut_sieie[std::make_pair(kTight,kEE)] = 0.031;

  cut_pfcharged[std::make_pair(kLoose,kEB)] = 2.6;
  cut_pfcharged[std::make_pair(kMedium,kEB)] = 1.5;
  cut_pfcharged[std::make_pair(kTight,kEB)] = 0.7;

  cut_pfcharged[std::make_pair(kLoose,kEE)] = 2.3;
  cut_pfcharged[std::make_pair(kMedium,kEE)] = 1.2;
  cut_pfcharged[std::make_pair(kTight,kEE)] = 0.5;

  cut_pfneutral[std::make_pair(kLoose,kEB)] = 3.5;
  cut_pfneutral[std::make_pair(kMedium,kEB)] = 1.0;
  cut_pfneutral[std::make_pair(kTight,kEB)] = 0.4;

  cut_pfneutral[std::make_pair(kLoose,kEE)] = 2.9;
  cut_pfneutral[std::make_pair(kMedium,kEE)] = 1.5;
  cut_pfneutral[std::make_pair(kTight,kEE)] = 1.5;

  cut_pfphoton[std::make_pair(kLoose,kEB)] = 1.3;
  cut_pfphoton[std::make_pair(kMedium,kEB)] = 0.7;
  cut_pfphoton[std::make_pair(kTight,kEB)] = 0.5;

  cut_pfphoton[std::make_pair(kLoose,kEE)] = 9999;
  cut_pfphoton[std::make_pair(kMedium,kEE)] = 1.0;
  cut_pfphoton[std::make_pair(kTight,kEE)] = 1.0;
  
}

bool StandardPhotonID::passID(const VecbosPho& pho,WP wp) {
  Region reg = (fabs(pho.SC.eta) <1.48 ? kEB:kEE);

  if(pho.HoverE > cut_HE[std::make_pair(wp,reg)]) return false;
  if(pho.SC.sigmaIEtaIEta > cut_sieie[std::make_pair(wp,reg)]) return false;
  
  return true;
}

bool StandardPhotonID::passIso(const VecbosPho& pho, WP wp) {
  Region reg = (fabs(pho.SC.eta) <1.48 ? kEB:kEE);
  float pt = pho.energy/cosh(pho.eta);

  if(pho.dr03ChargedHadronPFIso[0] > cut_pfcharged[std::make_pair(wp,reg)]) return false;
  if(pho.dr03NeutralHadronPFIso > cut_pfneutral[std::make_pair(wp,reg)]+0.04*pt) return false;
  if(pho.dr03PhotonPFIso > cut_pfphoton[std::make_pair(wp,reg)]+0.005*pt) return false;

  return true;
}
