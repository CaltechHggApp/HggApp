#ifndef StandardPhotonID_hh
#define StandardPhotonID_hh

#include "VecbosEGObject.hh"

#include <map>
#include <bitset>

class StandardPhotonID{
public:
  StandardPhotonID();

  enum WP:unsigned int{kLoose,kMedium,kTight};
  enum Region:unsigned int{kEB,kEE};

  bool passID(const VecbosPho& pho,WP wp,float rho);
  bool passIso(const VecbosPho& pho,WP wp,float rho);
  std::bitset<5> cutResults(const VecbosPho& pho,WP wp,float rho);
  typedef std::map<std::pair<StandardPhotonID::WP,StandardPhotonID::Region>,float> cutMap;
protected:
  cutMap cut_HE;
  cutMap cut_sieie;
  cutMap cut_pfcharged;
  cutMap cut_pfneutral;
  cutMap cut_pfphoton;

  std::vector<std::pair<float,float>> chEA;
  std::vector<std::pair<float,float>> nhEA;
  std::vector<std::pair<float,float>> phEA;
};

#endif
