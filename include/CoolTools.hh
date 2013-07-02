
//
// Written by C. Rogan
// Caltech
// 23-06-08
//
// Big ups to Numerical Recipes in C (minim, jacobi transform etc.)
//

#ifndef CoolTools_h
#define CoolTools_h

// std includes
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include <numeric>
#include <list>
#include "combination.hh"

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TVector.h>
#include <TLorentzVector.h>
#include <TMath.h>

// VecbosApp includes
#include "Jet.hh"
#include "MET.hh"
#include "CaloTower.hh"
#include "VecbosBase.hh"

using namespace std;

class CoolTools {

public:
  CoolTools();
  virtual ~CoolTools();

  vector<Jet> BoostJets(vector<Jet> input, TVector3 b);

  vector<TLorentzVector> BoostVectors(vector<TLorentzVector> input,
                                      TVector3 b);

  vector<Jet> CaloTowers2Jets(vector<CaloTower>, int);

  vector<TLorentzVector> Get4Vectors(vector<Jet> input);

  void CalcSphericity(vector<TLorentzVector> input);

  void CalcTSphericity(vector<TLorentzVector> input);

  void CalcThrust(vector<TLorentzVector> input);

  void CalcTranThrust(vector<TLorentzVector> input);

  double Legendre(double x, int l, int m);

  double FoxWolfram(vector<TLorentzVector>, int l);

  double TranFoxWolfram(vector<TLorentzVector>, int l);
  
  double FoxWolfram_test(vector<Jet>, int l);
  
  void SphericalHarmonic(double cos_theta, double phi, int l, int m,
			 double& real, double& imag);

  double Aplanarity() {return (1.5*S_e3);}
  
  double Sphericity() {return (1.5*(S_e2+S_e3));}

  TVector3 SphericityAxis() {return S_v1;}

  double detTS() {return (ST_e1*ST_e2);}

  double TranSphericity() {return (2.0*ST_e2/(ST_e1+ST_e2));}

  TVector3 TSphericityAxis() {return ST_v1;}

  double Thrust() {return T_1;}

  TVector3 ThrustAxis() {return T_v1;}

  TVector3 TranThrustAxis() {return TT_v;}

  TVector3 ThrustMajorAxis() {return T_v2;}

  TVector3 ThrustMinorAxis() {return T_v3;}

  double Oblateness() {return (T_2-T_3);}
  
  double ThrustMajor() {return T_2;}
  
  double ThrustMinor() {return T_3;}

  double TranThrust() {return TT;}
  
  double TranThrustMinor() {return TTm;}

  /// 3D CLEO cones w.r.t the Sphericity Axis
  void MakeCLEOCones_SphAxis(int ncons, vector<TLorentzVector> input);
  /// 3D CLEO cones w.r.t the Thrust Axis
  void MakeCLEOCones_ThrAxis(int ncons, vector<TLorentzVector> input);
  /// Transverse CLEO cones w.r.t the transverse Sphericity Axis
  void MakeTCLEOCones_SphAxis(int ncons, vector<TLorentzVector> input);
  /// Transverse CLEO cones w.r.t the transverse Thrust Axis
  void MakeTCLEOCones_ThrAxis(int ncons, vector<TLorentzVector> input);
  /// Transverse CLEO cones w.r.t the pT component of a generic axis
  void MakeTCLEOCones(int ncons, vector<TLorentzVector> input, TVector3 axis);
  /// 3D CLEO cones w.r.t a generic axis
  void MakeCLEOCones(int ncons, vector<TLorentzVector> input, TVector3 axis);
  /// Get The momentum in CLEO cone # ncon
  double GetMomCLEOCones(int ncon);
  /// Get the total momentum flow
  double GetMomFlowCLEOCones();


private:
  double S_e1, S_e2, S_e3;
  TVector3 S_v1, S_v2, S_v3;

  double ST_e1, ST_e2;
  TVector3 ST_v1, ST_v2;

  double T_1, T_2, T_3;
  TVector3 T_v1, T_v2, T_v3; // thrust axis, major and minor resp.

  double TT, TTm;
  TVector3 TT_v;

  TVector3 temp_rot_1;
  TVector3 temp_rot_2;

  //Calculation of Thrust globals
  double T_pcom[2], T_xicom[2];
  vector<TLorentzVector> T_list;
  
  //Thrust related functions
  double T_df1dim(double x);
  double T_f1dim(double x);
  void T_mnbrak(double *ax, double *bx, double *cx, 
		double *fa, double *fb, double *fc);
  double T_func(double p[2]);
  void T_dfunc(double x[2], double dfx[2]);
  double T_dlinmin(double p[2], double xi[2]);
  double T_dbrent(double ax, double bx, double cx, double tol, double *xmin);
  double TT_brent(double tol);
  double TT_func(double px);
  void TT_mnbrak(double *ax, double *bx, double *cx, 
		 double *fa, double *fb, double *fc);
  //for Sphericity calculation
  void jacobi(double a[3][3], int n, double d[3], double v[3][3]);

  // CLEO CONES
  vector<double>_eflowarray;
  double _momsum;
};




#endif
