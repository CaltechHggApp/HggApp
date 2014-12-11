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


  //effective areas
  chEA.push_back(std::make_pair(1.0,0.012));
  chEA.push_back(std::make_pair(1.479,0.010));
  chEA.push_back(std::make_pair(2.0,0.014));
  chEA.push_back(std::make_pair(2.2,0.012));
  chEA.push_back(std::make_pair(2.3,0.016));
  chEA.push_back(std::make_pair(2.4,0.020));
  chEA.push_back(std::make_pair(999,0.012));

  nhEA.push_back(std::make_pair(1.0,0.030));
  nhEA.push_back(std::make_pair(1.479,0.057));
  nhEA.push_back(std::make_pair(2.0,0.039));
  nhEA.push_back(std::make_pair(2.2,0.015));
  nhEA.push_back(std::make_pair(2.3,0.024));
  nhEA.push_back(std::make_pair(2.4,0.039));
  nhEA.push_back(std::make_pair(999,0.072));

  phEA.push_back(std::make_pair(1.0,0.148));
  phEA.push_back(std::make_pair(1.479,0.130));
  phEA.push_back(std::make_pair(2.0,0.112));
  phEA.push_back(std::make_pair(2.2,0.216));
  phEA.push_back(std::make_pair(2.3,0.262));
  phEA.push_back(std::make_pair(2.4,0.260));
  phEA.push_back(std::make_pair(999,0.266));

}

bool StandardPhotonID::passID(const VecbosPho& pho,WP wp,float rho) {
  std::bitset<5> res = cutResults(pho,wp,rho);
  if(res[0] || res[1]) return false; // two least significant bits are ID
  return true;
}

bool StandardPhotonID::passIso(const VecbosPho& pho, WP wp,float rho) {
  std::bitset<5> res = cutResults(pho,wp,rho);
  if(res[2] || res[3] || res[4]) return false; // three most significant bits are iso
  return true;

}

std::bitset<5> StandardPhotonID::cutResults(const VecbosPho& pho, WP wp,float rho) {
  Region reg = (fabs(pho.SC.eta) <1.48 ? kEB:kEE);
  float pt = pho.energy/cosh(pho.eta);

  float eta = fabs(pho.SC.eta);

  float thisCHEA = 0;
  for(auto ea: chEA) {
    if(ea.first > eta) {
      thisCHEA = ea.second;
      break;
    }
  }

  float thisNHEA = 0;
  for(auto ea: nhEA) {
    if(ea.first > eta) {
      thisNHEA = ea.second;
      break;
    }
  }

  float thisPHEA = 0;
  for(auto ea: phEA) {
    if(ea.first > eta) {
      thisPHEA = ea.second;
      break;
    }
  }

  std::bitset<5> res;
  if(pho.HTowOverE > cut_HE[std::make_pair(wp,reg)]) res.set(0);
  if(pho.SC.sigmaIEtaIEta > cut_sieie[std::make_pair(wp,reg)]) res.set(1);
  if(pho.dr03ChargedHadronPFIso[0] - rho*thisCHEA > cut_pfcharged[std::make_pair(wp,reg)]) res.set(2);
  if(pho.dr03NeutralHadronPFIso - rho*thisNHEA > cut_pfneutral[std::make_pair(wp,reg)]+0.04*pt) res.set(3);
  if(pho.dr03PhotonPFIso - rho*thisPHEA > cut_pfphoton[std::make_pair(wp,reg)]+0.005*pt) res.set(4);
  
  return res;
}
