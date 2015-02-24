#include "TRandom3.h"


#include <vector>
class evalError {
public:
  evalError() : rng(0) {}
  
  void addError(float e) { errors.push_back(e); }
  void addPercentError(float e) { percentErrors.push_back(e); }
  
  float eval(float val);
  void clear() { errors.clear(); percentErrors.clear(); }

protected:
  TRandom3 rng;

  std::vector<float> errors;
  std::vector<float> percentErrors;
};

float evalError::eval(float val) {
  float offset=0;

  for( int iErr=0; iErr<errors.size(); iErr++) {
    float frac = errors.at(iErr)/val;
    float thisOff = rng.Gaus(0,frac);
    offset+=thisOff;
  }
  for( int iErr=0; iErr<percentErrors.size(); iErr++) {
    float frac = percentErrors.at(iErr);
    float thisOff = rng.Gaus(0,frac);
    offset+=thisOff;
  }
  
  return val*(1+offset);
}
