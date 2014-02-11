#include "StandardElectronID.hh"


StandardElectronID::StandardElectronID() {
  cut_dEtaIn[std::make_pair(kVeto,kEB)] = 0.007;
  cut_dEtaIn[std::make_pair(kLoose,kEB)] = 0.007;
  cut_dEtaIn[std::make_pair(kMedium,kEB)] = 0.004;
  cut_dEtaIn[std::make_pair(kTight,kEB)] = 0.004;

  cut_dPhiIn[std::make_pair(kVeto,kEB)] = 0.8;
  cut_dPhiIn[std::make_pair(kLoose,kEB)] = 0.15;
  cut_dPhiIn[std::make_pair(kMedium,kEB)] = 0.06;
  cut_dPhiIn[std::make_pair(kTight,kEB)] = 0.03;

  cut_sieie[std::make_pair(kVeto,kEB)] = 0.01;
  cut_sieie[std::make_pair(kLoose,kEB)] = 0.01;
  cut_sieie[std::make_pair(kMedium,kEB)] = 0.01;
  cut_sieie[std::make_pair(kTight,kEB)] = 0.01;

  cut_HE[std::make_pair(kVeto,kEB)] = 0.15;
  cut_HE[std::make_pair(kLoose,kEB)] = 0.12;
  cut_HE[std::make_pair(kMedium,kEB)] = 0.12;
  cut_HE[std::make_pair(kTight,kEB)] = 0.12;

  cut_d0Track[std::make_pair(kVeto,kEB)] = 0.04;
  cut_d0Track[std::make_pair(kLoose,kEB)] = 0.02;
  cut_d0Track[std::make_pair(kMedium,kEB)] = 0.02;
  cut_d0Track[std::make_pair(kTight,kEB)] = 0.02;

  cut_dzTrack[std::make_pair(kVeto,kEB)] = 0.2;
  cut_dzTrack[std::make_pair(kLoose,kEB)] = 0.2;
  cut_dzTrack[std::make_pair(kMedium,kEB)] = 0.1;
  cut_dzTrack[std::make_pair(kTight,kEB)] = 0.1;

  cut_1oE1oP[std::make_pair(kVeto,kEB)] = 9999;
  cut_1oE1oP[std::make_pair(kLoose,kEB)] = 0.05;
  cut_1oE1oP[std::make_pair(kMedium,kEB)] = 0.05;
  cut_1oE1oP[std::make_pair(kTight,kEB)] = 0.05;

  cut_PFIso[std::make_pair(kVeto,kEB)] = 0.15;
  cut_PFIso[std::make_pair(kLoose,kEB)] = 0.15;
  cut_PFIso[std::make_pair(kMedium,kEB)] = 0.15;
  cut_PFIso[std::make_pair(kTight,kEB)] = 0.10;

  cut_missingHit[std::make_pair(kVeto,kEB)] = 9999;
  cut_missingHit[std::make_pair(kLoose,kEB)] = 1.5; //these are really integers, so tlets not have any float->int conversion issues
  cut_missingHit[std::make_pair(kMedium,kEB)] = 1.5;
  cut_missingHit[std::make_pair(kTight,kEB)] = 0.5;

  //EE

  cut_dEtaIn[std::make_pair(kVeto,kEE)] = 0.01;
  cut_dEtaIn[std::make_pair(kLoose,kEE)] = 0.009;
  cut_dEtaIn[std::make_pair(kMedium,kEE)] = 0.007;
  cut_dEtaIn[std::make_pair(kTight,kEE)] = 0.005;

  cut_dPhiIn[std::make_pair(kVeto,kEE)] = 0.7;
  cut_dPhiIn[std::make_pair(kLoose,kEE)] = 0.10;
  cut_dPhiIn[std::make_pair(kMedium,kEE)] = 0.03;
  cut_dPhiIn[std::make_pair(kTight,kEE)] = 0.02;

  cut_sieie[std::make_pair(kVeto,kEE)] = 0.03;
  cut_sieie[std::make_pair(kLoose,kEE)] = 0.03;
  cut_sieie[std::make_pair(kMedium,kEE)] = 0.03;
  cut_sieie[std::make_pair(kTight,kEE)] = 0.03;

  cut_HE[std::make_pair(kVeto,kEE)] = 99999;
  cut_HE[std::make_pair(kLoose,kEE)] = 0.10;
  cut_HE[std::make_pair(kMedium,kEE)] = 0.10;
  cut_HE[std::make_pair(kTight,kEE)] = 0.10;

  cut_d0Track[std::make_pair(kVeto,kEE)] = 0.04;
  cut_d0Track[std::make_pair(kLoose,kEE)] = 0.02;
  cut_d0Track[std::make_pair(kMedium,kEE)] = 0.02;
  cut_d0Track[std::make_pair(kTight,kEE)] = 0.02;

  cut_dzTrack[std::make_pair(kVeto,kEE)] = 0.2;
  cut_dzTrack[std::make_pair(kLoose,kEE)] = 0.2;
  cut_dzTrack[std::make_pair(kMedium,kEE)] = 0.1;
  cut_dzTrack[std::make_pair(kTight,kEE)] = 0.1;

  cut_1oE1oP[std::make_pair(kVeto,kEE)] = 9999;
  cut_1oE1oP[std::make_pair(kLoose,kEE)] = 0.05;
  cut_1oE1oP[std::make_pair(kMedium,kEE)] = 0.05;
  cut_1oE1oP[std::make_pair(kTight,kEE)] = 0.05;

  cut_PFIso[std::make_pair(kVeto,kEE)] = 0.15;
  cut_PFIso[std::make_pair(kLoose,kEE)] = 0.15;
  cut_PFIso[std::make_pair(kMedium,kEE)] = 0.15;
  cut_PFIso[std::make_pair(kTight,kEE)] = 0.10;

  cut_missingHit[std::make_pair(kVeto,kEE)] = 9999;
  cut_missingHit[std::make_pair(kLoose,kEE)] = 1.5; //these are really integers, so tlets not have any float->int conversion issues
  cut_missingHit[std::make_pair(kMedium,kEE)] = 1.5;
  cut_missingHit[std::make_pair(kTight,kEE)] = 0.5;
}

bool StandardElectronID::passID(const VecbosEle& ele, WP wp) {
  std::bitset<__NUM_ELE_ID_CUTS__> res = cutResults(ele,wp);
  
  return ! (res&std::bitset<__NUM_ELE_ID_CUTS__>(0xFF)).any(); // 8 lowest bits are the ID bits
}

bool StandardElectronID::passIso(const VecbosEle& ele, WP wp) {
  std::bitset<__NUM_ELE_ID_CUTS__> res = cutResults(ele,wp);
  
  return ! (res&std::bitset<__NUM_ELE_ID_CUTS__>(0x100)).any(); // highest bit is the Iso Bit
}

bool StandardElectronID::passAll(const VecbosEle& ele, WP wp) {
  std::bitset<__NUM_ELE_ID_CUTS__> res = cutResults(ele,wp);
  
  return ! (res.any());
}




std::bitset<__NUM_ELE_ID_CUTS__> StandardElectronID::cutResults(const VecbosEle& ele, WP wp) {
  std::bitset<__NUM_ELE_ID_CUTS__> res;
  Region r = kEB;
  if(fabs(ele.SC.eta)>1.48) r = kEE;

  if(ele.dEtaSCTrackAtVtx > cut_dEtaIn[std::make_pair(wp,r)]) res.set(0);
  if(ele.dPhiSCTrackAtVtx > cut_dPhiIn[std::make_pair(wp,r)]) res.set(1);
  if(ele.SC.sigmaIEtaIEta > cut_sieie[std::make_pair(wp,r)]) res.set(2);
  if(ele.HoverE > cut_HE[std::make_pair(wp,r)]) res.set(3);
  if(ele.d0Track > cut_d0Track[std::make_pair(wp,r)]) res.set(4);
  if(ele.dzTrack > cut_dzTrack[std::make_pair(wp,r)]) res.set(5);
  float eMp = 1./ele.SC.energy*fabs(1. - ele.EOverP);
  if(eMp > cut_1oE1oP[std::make_pair(wp,r)]) res.set(6);
  if(ele.expInnerLayersHits > cut_missingHit[std::make_pair(wp,r)]) res.set(7);

  float pt = ele.pt;
  float pfiso = ele.dr03ChargedHadronPFIso+ele.dr03NeutralHadronPFIso+ele.dr03PhotonPFIso;
  pfiso/=pt;
  if(pfiso > cut_PFIso[std::make_pair(wp,r)]) res.set(8);

  return res;
}
