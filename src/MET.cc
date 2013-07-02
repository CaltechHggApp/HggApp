// VecbosApp includes
#include "include/MET.hh"
  
MET::MET(double mex, double mey, double mez, double sumet){
  et = sqrt(mex*mex + mey*mey);
  ex = mex;
  ey = mey;
  Phi = atan2(mex, mey);
  Sumet = sumet;
  ez = mez;
}

MET::~MET() {}
