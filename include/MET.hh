#ifndef MET_h
#define MET_h

// ROOT includes
#include <TLorentzVector.h>

/// MET class, to represent missing transverse energy.
/// Methods are self-documenting.
class MET {
public:
  /// Class Constructor
  MET(double mex, double mey, double mez, double sumet);
  /// Class Destructor
  virtual ~MET();

  double met() {return et;}  
  double mex() {return ex;}  
  double mey() {return ey;}  
  double mez() {return ez;}  
  double phi() {return Phi;} 
  double sumet() {return Sumet;} //< Returns ... 
  TLorentzVector metVector() {TLorentzVector v(ex, ey, 0.0, et); return v;} //< Returns TLorentzVector of MET.

private:
  double et;
  double ex;
  double ey;
  double ez;
  double Phi;
  double Sumet;
};

#endif
