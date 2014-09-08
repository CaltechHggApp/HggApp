#ifndef RazorVariables_hh
#define RazorVariables_hh

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "exception"

class TooManyJets : public std::exception {
public:
  virtual  const char* what() const throw() {
    return "Trying to combine too many jets into hemispheres, this is a FATAL error";
  }
};

class TooFewJets : public std::exception {
public: 
  virtual const char* what() const throw() {
    return "Trying to combine too few jets into hemispheres.  This error should be caught and handled";
  }
};

class RazorVariables {
public:
  RazorVariables(){}

  static void CombineJets(const std::vector<TLorentzVector>& jets, TLorentzVector& hem1, TLorentzVector& hem2,std::vector<int>* hemAssignment=0) throw(TooManyJets,TooFewJets);
  static double CalcGammaMRstar(const TLorentzVector& hem1, const TLorentzVector& hem2) throw();
  static double CalcMTR(const TLorentzVector&hem1, const TLorentzVector&hem2,const TVector3& met) throw();
};

#endif
