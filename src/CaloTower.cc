// VecbosApp includes
#include "include/CaloTower.hh"

CaloTower::CaloTower(double Em, double Had, TVector3 calo, TVector3 ecal,TVector3 hcal){
  double E0 = Em + Had;
  if(E0 > 0.0){
    Efrac = Em/E0;
    Hfrac = Had/E0;
  } else {
    Efrac = 0.0;
    Hfrac = 0.0;
  }
  ECAL = ecal;
  CALO = calo;
  HCAL = hcal;
  TVector3 v2;
  
  v2.SetPtEtaPhi((1.0/cosh(calo.Eta()))*E0, calo.Eta(), calo.Phi());
 
  v.SetPtEtaPhiE(v2.Pt(), v2.Eta(), v2.Phi(), v2.Mag());
  
}


CaloTower::~CaloTower() {}

