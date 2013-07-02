#define _onePi 3.141592653
#define _twoPi 6.283185307
#include <TChain.h>
#include <math.h>
////using namespace math;

struct ECAL_GEO{
  // EB gap positions
  double _barrelCGap[169][360][2];
  double _barrelSGap[33][180][2];
  double _barrelTGap[33][180][2];
  double _barrelMGap[7][18][2];

  // EE crystal existence and gap positions
  bool   _endcapCrystal[100][100];
  double _endcapCGap[2][7080][2];
  double _endcapSGap[2][264][2];
  double _endcapMGap[2][1][2];

  double _endcapCGapxy[2][7080][2];
  double _endcapSGapxy[2][264][2];
  double _endcapMGapxy[2][1][2];
};

struct THIS_ECAL_GEO{
  // Actual data for each instantiated object
  unsigned _be,_hl;
  double _e,_eta,_phi,_r9;
  double _aC,_aS,_aM,_bC,_bS,_bM;
  double _aT,_bT;
  double _xyaC,_xyaS,_xyaM,_xybC,_xybS,_xybM;
  double _xZ,_yZ; 
};

ECAL_GEO loadecalGapCoordinates(Bool_t isRealData){
  //load the ECAL geometry (gap coordinates) from a file and return it as an object
  TChain *feg = new TChain("ecalGap");

  if(isRealData){
    feg->Add("/home/amott/HggApp/ecalGap.cmssw425.GeometryDB.FT_R_42_V17A.root");///data 
  }else{
    ///feg->Add("/afs/cern.ch/user/y/yangyong/w0/ecalGap.cmssw425.root"); //mc
    feg->Add("/home/amott/HggApp/ecalGap.cmssw425.GeometryDB.START42_V12_new.root"); //MC
  }
  
  ECAL_GEO geometry;
  
  feg->SetBranchAddress("_barrelCGap", geometry._barrelCGap);
  feg->SetBranchAddress("_barrelSGap", geometry._barrelSGap);
  feg->SetBranchAddress("_barrelTGap", geometry._barrelTGap);
  feg->SetBranchAddress("_barrelMGap", geometry._barrelMGap);
  feg->SetBranchAddress("_endcapCGap", geometry._endcapCGap);
  feg->SetBranchAddress("_endcapSGap", geometry._endcapSGap);
  feg->SetBranchAddress("_endcapMGap", geometry._endcapMGap);
  feg->SetBranchAddress("_endcapCGapxy", geometry._endcapCGapxy);
  feg->SetBranchAddress("_endcapSGapxy", geometry._endcapSGapxy);
  feg->SetBranchAddress("_endcapMGapxy", geometry._endcapMGapxy);

  
  feg->GetEntry(0);
    
  return geometry;
}

double dPhi(double f0, double f1) {
  double df(f0-f1);
  if(df> _onePi) df-=_twoPi;
  if(df<-_onePi) df+=_twoPi;
  return df;
}

double aPhi(double f0, double f1) {
  double af(0.5*(f0+f1));
  if(fabs(dPhi(af,f0))>0.5*_onePi) {
    if(af>=0.0) af-=_onePi;
    else        af+=_onePi;
  }
  
  //  assert(fabs(dPhi(af,f0))<0.5*_onePi);
  //assert(fabs(dPhi(af,f1))<0.5*_onePi);
  
  return af;
}


double xZ(double eta, double phi){
  return asinh(cos(phi)/sinh(eta));
}

double yZ(double eta, double phi) {
  return asinh(sin(phi)/sinh(eta));
}



THIS_ECAL_GEO getGapCoordinates(ECAL_GEO geometry, double eta, double phi){
  // Check constants have been set up
  //assert(_initialised);
  
  THIS_ECAL_GEO thisGeometry;

  // Determine if EB or EE
  bool _be=(fabs(eta)<1.482?0:1);
  
  //  // Determine if high or low R9
  //   if(_be==0) _hl=(_r9>=0.94?0:1);
  //   else       _hl=(_r9>=0.95?0:1);
  
  // Coordinates relative to cracks
  double r2Min;
  if(_be==0) {
    
    r2Min=1.0e6;
    for(unsigned i(0);i<169;i++) {
      for(unsigned j(0);j<360;j++) {
	double de(eta-geometry._barrelCGap[i][j][0]);
	double df(dPhi(phi,geometry._barrelCGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=84) {
	    thisGeometry._aC= de;
	    thisGeometry._bC=-df;
	  } else {
	    thisGeometry._aC=-de;
	    thisGeometry._bC= df;
	  }
	}
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<33;i++) {
      for(unsigned j(0);j<180;j++) {
	double de(eta-geometry._barrelSGap[i][j][0]);
	double df(dPhi(phi,geometry._barrelSGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=16) {
	    thisGeometry._aS= de;
	    thisGeometry._bS=-df;
	    
	  } else {
	    thisGeometry._aS=-de;
	    thisGeometry._bS= df;
	  }
	}
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<7;i++) {
      for(unsigned j(0);j<18;j++) {
	double de(eta-geometry._barrelMGap[i][j][0]);
	double df(dPhi(phi,geometry._barrelMGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=3) {
	    thisGeometry._aM= de;
	    thisGeometry._bM=-df;
	  } else {
	    thisGeometry._aM=-de;
	    thisGeometry._bM= df;
	  }
	}
      }
    }
    
    
    r2Min=1.0e6;
    for(unsigned i(0);i<33;i++) {
      for(unsigned j(0);j<72;j++) {
	double de(eta-geometry._barrelTGap[i][j][0]);
	double df(dPhi(phi,geometry._barrelTGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=16) {
	    thisGeometry._aT= de;
	    thisGeometry._bT=-df;
	  } else {
	    thisGeometry._aT=-de;
	    thisGeometry._bT= df;
	  }
	}
      }
    }
    
  } else {
    unsigned iz(eta>=0.0?0:1);

    //double r[2]={xZ(eta,phi),yZ(eta,phi)};
    
    thisGeometry._xZ = xZ(eta,phi);
    thisGeometry._yZ = yZ(eta,phi);
    double r[2] = {thisGeometry._xZ, thisGeometry._yZ};
    
    r2Min=1.0e6;
    for(unsigned i(0);i<7080;i++) {
      double dx(r[0]-geometry._endcapCGap[iz][i][0]);
      double dy(r[1]-geometry._endcapCGap[iz][i][1]);
      double r2(dx*dx+dy*dy);

      if(r2<r2Min) {
	r2Min=r2;
	if(r[0]>0.0) thisGeometry._aC= dx;
	else         thisGeometry._aC=-dx;
	if(r[1]>0.0) thisGeometry._bC= dy;
	else         thisGeometry._bC=-dy;
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<264;i++) {
      double dx(r[0]-geometry._endcapSGap[iz][i][0]);
      double dy(r[1]-geometry._endcapSGap[iz][i][1]);
      double r2(dx*dx+dy*dy);

      if(r2<r2Min) {
	r2Min=r2;
	if(r[0]>0.0) thisGeometry._aS= dx;
	else         thisGeometry._aS=-dx;
	if(r[1]>0.0) thisGeometry._bS= dy;
	else         thisGeometry._bS=-dy;
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<1;i++) {
      double dx(r[0]-geometry._endcapMGap[iz][i][0]);
      double dy(r[1]-geometry._endcapMGap[iz][i][1]);
      double r2(dx*dx+dy*dy);

      if(r2<r2Min) {
	r2Min=r2;
	if(iz==0) {thisGeometry._aM= dx;thisGeometry._bM= dy;}
	else      {thisGeometry._aM=-dx;thisGeometry._bM=-dy;}
      }
    }
  }

  return thisGeometry;
}



THIS_ECAL_GEO getGapCoordinatesXY(ECAL_GEO geometry, double x, double y,double eta){
  // Check constants have been set up
  //assert(_initialised);
  
  THIS_ECAL_GEO thisGeometry;
  // Determine if EB or EE
  ///_be=(fabs(eta)<1.482?0:1);
  
  if( fabs(eta) < 1.482) return thisGeometry;  //for endcap only
  
  unsigned iz(eta>=0.0?0:1);
  double r[2]={x,y};
  
  double r2Min=1.0e6;
  for(unsigned i(0);i<7080;i++) {
    double dx(r[0]-geometry._endcapCGapxy[iz][i][0]);
    double dy(r[1]-geometry._endcapCGapxy[iz][i][1]);
    double r2(dx*dx+dy*dy);

    if(r2<r2Min) {
      r2Min=r2;
      if(r[0]>0.0) thisGeometry._xyaC= dx;
      else         thisGeometry._xyaC=-dx;
      if(r[1]>0.0) thisGeometry._xybC= dy;
      else         thisGeometry._xybC=-dy;
    }
  }
    
  r2Min=1.0e6;
  for(unsigned i(0);i<264;i++) {
    double dx(r[0]-geometry._endcapSGapxy[iz][i][0]);
    double dy(r[1]-geometry._endcapSGapxy[iz][i][1]);
    double r2(dx*dx+dy*dy);

    if(r2<r2Min) {
      r2Min=r2;
      if(r[0]>0.0) thisGeometry._xyaS= dx;
      else         thisGeometry._xyaS=-dx;
      if(r[1]>0.0) thisGeometry._xybS= dy;
      else         thisGeometry._xybS=-dy;
    }
  }
    
  r2Min=1.0e6;
  for(unsigned i(0);i<1;i++) {
    double dx(r[0]-geometry._endcapMGapxy[iz][i][0]);
    double dy(r[1]-geometry._endcapMGapxy[iz][i][1]);
    double r2(dx*dx+dy*dy);

    if(r2<r2Min) {
      r2Min=r2;
      if(iz==0) {thisGeometry._xyaM= dx;thisGeometry._xybM= dy;}
      else      {thisGeometry._xyaM=-dx;thisGeometry._xybM=-dy;}
    }
  }

  return thisGeometry;
}

