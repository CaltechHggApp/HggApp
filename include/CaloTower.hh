#ifndef CaloTower_h
#define CaloTower_h

// ROOT includes
#include <TVector3.h>
#include <TLorentzVector.h>

/// CaloTower class, to represent a CaloTower.
/// Methods are self-documenting.
class CaloTower {
public:

  /// Class Constructor
  CaloTower(double Em, double Had, TVector3 calo, TVector3 ecal, TVector3 hcal); 
  /// Class Destructor
  virtual ~CaloTower();     

  double px(){return v.Px();}
  double py(){return v.Py();}
  double pz(){return v.Pz();}
  double e(){return v.E();}

  double et() {return v.Et();}
  double eta() {return v.Eta();}
  double phi() {return v.Phi();}
  double theta() {return v.Theta();}
  double EmFrac() {return Efrac;} //< Returns the eletromagnetic fraction of the calotower.
  double HadFrac() {return Hfrac;} //< Returns the hadronic fraction of the calotower.
  double EmEnergy() {double val = Efrac*v.E(); return val;} //< Returns EmFrac()*e()
  double HadEnergy() {double val = Hfrac*v.E(); return val;} //<Returns HadFrac*e()
  TVector3 getCALOPos() {return CALO;} //< Returns the 3D position (TVector3) of the calotower. 
  TVector3 getECALPos() {return ECAL;} //< Returns the 3D position (TVector3) of the ECAL part of the calotower.
  TVector3 getHCALPos() {return HCAL;} //< Returns the 3D position (TVector3) of the HCAL part of the calotower.
  TLorentzVector Get4Vector() {return v;}
  TVector3 Get3Vector() {return v.Vect();}

private:
  TLorentzVector v;
  TVector3 ECAL;
  TVector3 HCAL;
  TVector3 CALO;
  double Efrac;
  double Hfrac;
};

#endif
