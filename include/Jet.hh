#ifndef Jet_h
#define Jet_h

// ROOT includes
#include <TLorentzVector.h>

/// Jet class, to represent a jet.
/// Methods are self-documenting.
class Jet {
public:
  /// Class Constructor
  Jet(TLorentzVector p, double efrac, double hfrac);
  /// Class Destructor
  virtual ~Jet();

  void setPFJetID(double nHFrac, double nEmFrac, int nConst, double chHFrac, double chMult, double chEmFrac);

  double mass() {return v.M();}
  double et() {return v.Et();} 
  double phi() {return v.Phi();}
  double eta() {return v.Eta();}
  double theta() {return v.Theta();}
  double pt() {return v.Pt();}
  double px() {return v.Px();}
  double py() {return v.Py();}
  double pz() {return v.Pz();}
  double e() {return v.E();}
  double p() {return v.P();}

  /// Thiago: Another set of functions to define a common interface with TLorentzVector.
  /// Given this, perhaphs now we should make Jet inherit from TLorentzVector?
  double M() {return v.M();}
  double Et() {return v.Et();} 
  double Phi() {return v.Phi();}
  double Eta() {return v.Eta();}
  double Theta() {return v.Theta();}
  double Pt() {return v.Pt();}
  double Px() {return v.Px();}
  double Py() {return v.Py();}
  double Pz() {return v.Pz();}
  double E() {return v.E();}
  double P() {return v.P();}

  double EmFrac() {return Efrac;} //< Returns the eletromagnetic fraction of the jet.
  double HadFrac() {return Hfrac;} //< Returns the hadronic fraction of the jet.

  TLorentzVector Sum(Jet other) {return v+other.v;}

  TLorentzVector Get4Vector() {return v;}
  TVector3 Get3Vector() {return v.Vect();}

  bool isPFJetID(int WP);

  enum jetIdWP { none=0, loose=1, medium=2, tight=3 };

private:
  TLorentzVector v;
  double Efrac;
  double Hfrac;

  // PFJet ID
  double neutralHadFrac, neutralEmFraction;
  int nConstituents;
  double chargedHadFraction, chargedMultiplicity, chargedEmFraction;

};

#endif
