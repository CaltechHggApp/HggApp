// std includes
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>

using namespace std;

// vecbosApp includes
#include "VertexProjector.hh"

VertexProjector::VertexProjector(double x0, double y0, double z0) {
  VtxX0 = x0;
  VtxY0 = y0;
  VtxZ0 = z0;
  EcalRad = 1.5; // [m], to check
}

VertexProjector::~VertexProjector() {}


void VertexProjector::SetLorentzVector(double thisE, double thispx, double thispy, double thispz) {
  E  = thisE;
  px = thispx;  
  py = thispy;
  pz = thispz;
  CalcNewVector();
}

void VertexProjector::CalcNewVector(){

  Theta = atan2(sqrt(px*px+py*py),pz);
  Phi = atan2(py,px);
  
  // Calculates the new direction
  // from the starting Direction R
  // and the new vertex Vtx = (VtxX0VtxY0,VtxZ0)
  // as R'=R-Vtx
  
  double Rx = EcalRad*cos(Phi) - VtxX0;
  double Ry = EcalRad*sin(Phi) - VtxY0;
  double Rz = EcalRad/tan(Theta) - VtxZ0;
  
  // Calculate the new angles
  Theta = atan2(sqrt(Rx*Rx+Ry*Ry),Rz);
  Phi   = atan2(Ry,Rx);
  
  // Calculate the new Lorents momentum
  // E is untached, since it is the measurement
  // We do not change p (which is strictly exact only
  // if p=E, as in the case of CaloTowers as input)
  double p = sqrt(px*px+py*py+pz*pz);

  px = p*sin(Theta)*cos(Phi);
  py = p*sin(Theta)*sin(Phi);
  pz = p*cos(Theta);

}

