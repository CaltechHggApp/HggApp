#ifndef StandardElectronID_hh
#define StandardElectronID_hh

#include "VecbosEGObject.hh"

#include <map>
#include <bitset>

#define __NUM_ELE_ID_CUTS__ 9

class StandardElectronID {
public:
  StandardElectronID();
  enum WP:unsigned int{kVeto=0,kLoose,kMedium,kTight};
  enum Region{kEB,kEE};
  
  bool passID(const VecbosEle& ele,WP wp);
  bool passIso(const VecbosEle& ele,WP wp);
  bool passAll(const VecbosEle& ele,WP wp);

  std::bitset<__NUM_ELE_ID_CUTS__> cutResults(const VecbosEle& ele,WP wp);
  typedef std::map<std::pair<StandardElectronID::WP,StandardElectronID::Region>,float> cutMap;
protected:
  cutMap cut_dEtaIn;
  cutMap cut_dPhiIn;
  cutMap cut_sieie;
  cutMap cut_HE;
  cutMap cut_d0Track;
  cutMap cut_dzTrack;
  cutMap cut_1oE1oP;
  cutMap cut_PFIso;
  cutMap cut_missingHit;
};

#endif
