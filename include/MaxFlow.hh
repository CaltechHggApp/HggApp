#ifndef _maxflow_hh_included_
#define _maxflow_hh_included_
#include <vector>
#include <TVector2.h>
#include <TFitter.h>

// Function to calculate numerator.
double calculateNumer(TVector2& thisVector, std::vector<TVector2>* v) {
  double res = 0.;
  int vsize = v->size();
  
  for(int i=0; i!=vsize; ++i)
    res += fabs(thisVector*v->at(i));
  
  // Returns -res because we are minimizing instead of maximizing.
  return (-res);
};

// Auxiliary function.
double myFunction(double phi) {
  extern std::vector<TVector2>* TracksPtr; 
  TVector2 tmp;
  tmp.SetMagPhi(1.0,phi);
  return calculateNumer(tmp, TracksPtr);
};

// Minuit function.
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  result = myFunction(par[0]);
};

#endif
