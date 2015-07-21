#include "TMinuit.h"
#include "TMath.h"

#include <iostream>

#include "TH1D.h"

//the parameters of the bin we are in
struct par {
  double obs=0;
  double Nupper=3,Nlower=0;
  double UpperSF=0.383,UpperSFe=0.0021;
  double LowerSF=0.279,LowerSFe=0.0014;
  double Higgs=0.06,HiggsErr=0.012;
  double S=0;
} params;


double gaus(double x, double m, double s) {
  return 1/sqrt(TMath::Pi()*s)*TMath::Exp( -1*(x-m)*(x-m)/s/s);
}

double pois(double N, double alpha) {
  return TMath::Power(alpha,N)*TMath::Exp(-1*alpha);
}


//nuisance parameters
//float UpperSF,LowerSF,B,Higgs;

// npar should be 4 (free) parameters
// UpperSF,LowerSF,upperB,HiggsB
void negLikelihood(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t flag) {
  if(params.S+par[0]*par[2]+par[3] <0) { 
    f=0; //truncate so we can scan negative S
    return;
  }
  f = 1;
  f *= TMath::Gaus(par[0],params.UpperSF,params.UpperSFe); // upper SF compatibility
  f *= TMath::Gaus(par[1],params.LowerSF,params.LowerSFe); // lower SF compatibility
  f *= TMath::PoissonI(params.Nupper,par[2]);                    // upper sideband compatibility
  f *= TMath::PoissonI(params.Nlower,par[2]*par[0]/par[1]);      // lower sideband compatibility
  f *= TMath::Gaus(par[3],params.Higgs,params.HiggsErr);        // Higgs background compatibility
  f *= TMath::PoissonI(params.obs,params.S+par[0]*par[2]+par[3]);            // S compatibility
  f*=-1;
}

float getR(TMinuit& m);

//just return the -2delta LL , best S

std::pair<float,float> profileNoHistMinuit(float step=0.01) {
  //std::cout << params.obs << " " << params.Nupper << " " << params.Nlower << " " << params.UpperSF << " " << params.UpperSFe << " " << std::endl;
  //std::cout << params.LowerSF << " " << params.LowerSFe << " " << params.Higgs << " " << params.HiggsErr << " " << params.UpperSFe << " " << std::endl;
  

  const int nPar=4;

  TMinuit min(nPar);
  Int_t ierflg = 0;

  min.SetFCN(negLikelihood);

  double pStart[nPar]={params.UpperSF,params.LowerSF,params.Nupper,params.Higgs};
  double pStep [nPar]={0.0001,0.0001,0.0001,0.0001};
  

  Double_t arglist[10];
  Double_t arglistS0[10];

  arglist[0]=5000;
  arglist[1]=2.;

  arglistS0[0]=5000;
  arglistS0[1]=0.;

  min.SetPrintLevel(-1);

  float curMin=99999999999; // store the minimum value of the likelihood found
  float bestS=0;
  const int nReqGreater=int(1./step); // number of point greater than the curMin required to declare global min
  int nGreater=0;          // number of points currently found greater than curMin

  float likeS0 = 0; //store the likelihood at S=0

  //do the S=0 likelihood

  const int nMax=1E6; // never exceed this many calls

  params.S=0;
  min.mnparm(0,"UpperSF",pStart[0],pStep[0],0,0,ierflg);
  min.mnparm(1,"LowerSF",pStart[1],pStep[1],0,0,ierflg);
  min.mnparm(2,"UpperB",pStart[2],pStep[2],0,0,ierflg);
  min.mnparm(3,"HiggsB",pStart[3],pStep[3],0,0,ierflg);

  min.mnexcm("MIGRAD",arglistS0,2,ierflg);
  min.mnexcm("MIGRAD",arglist,2,ierflg);

  likeS0 = getR(min);
  if(TMath::IsNaN(likeS0)) {
    std::cout << "FATAL ERROR: S=0 likelihood is NaN" << std::endl;
    return std::make_pair(0.,0.);
  }
  curMin = likeS0;

  //scan S
  int nRun=0;

  while(true) {
    if(nRun++ > nMax) return std::make_pair(0,0);
    params.S+=step;
    min.mnparm(0,"UpperSF",pStart[0],pStep[0],0,0,ierflg);
    min.mnparm(1,"LowerSF",pStart[1],pStep[1],0,0,ierflg);
    min.mnparm(2,"UpperB",pStart[2],pStep[2],0,0,ierflg);
    min.mnparm(3,"HiggsB",pStart[3],pStep[3],0,0,ierflg);

    //std::cout << params.S << std::endl;
    min.mnexcm("MIGRAD",arglistS0,2,ierflg);
    min.mnexcm("MIGRAD",arglist,2,ierflg);
    if(ierflg!=0) {
      min.mnparm(0,"UpperSF",pStart[0]+0.002,pStep[0],0,0,ierflg);
      min.mnparm(1,"LowerSF",pStart[1]-0.002,pStep[1],0,0,ierflg);
      min.mnparm(2,"UpperB",pStart[2]+1,pStep[2],0,0,ierflg);
      min.mnparm(3,"HiggsB",pStart[3]+0.01,pStep[3],0,0,ierflg);
      
      
      min.mnexcm("MIGRAD",arglistS0,2,ierflg);
      min.mnexcm("MIGRAD",arglist,2,ierflg);
      if(ierflg!=0) {
	//std::cout << "ERROR: " << ierflg << std::endl;
	//std::cout << "S: " << params.S << std::endl;
	//min.SetPrintLevel(1);
	min.mnexcm("MIGRAD",arglist,2,ierflg);
	continue;
      }      
    }

    float thisL = getR(min);
    //std::cout << thisL << std::endl;
    if(TMath::IsNaN(thisL)) {
      std::cout << "NaN!  S=" << params.S << std::endl;
      continue;
    }

    //update minimum as needed
    if(thisL<curMin) {
      curMin=thisL;
      nGreater=0;
      bestS = params.S;
    } else {
      nGreater++;
    }

    if(nGreater>nReqGreater) break; //check nGreater
  }

  //now check the negative S

  params.S=0;
  step*=-1; //swap the sign of step
  nGreater=0;
  while(true) {
    if(nRun++ > nMax) return std::make_pair(0,0);
    params.S+=step;
    min.mnparm(0,"UpperSF",pStart[0],pStep[0],0,0,ierflg);
    min.mnparm(1,"LowerSF",pStart[1],pStep[1],0,0,ierflg);
    min.mnparm(2,"UpperB",pStart[2],pStep[2],0,0,ierflg);
    min.mnparm(3,"HiggsB",pStart[3],pStep[3],0,0,ierflg);

    //std::cout << params.S << std::endl;
    min.mnexcm("MIGRAD",arglistS0,2,ierflg);
    min.mnexcm("MIGRAD",arglist,2,ierflg);
    if(ierflg!=0) {
      min.mnparm(0,"UpperSF",pStart[0]+0.002,pStep[0],0,0,ierflg);
      min.mnparm(1,"LowerSF",pStart[1]-0.002,pStep[1],0,0,ierflg);
      min.mnparm(2,"UpperB",pStart[2]+1,pStep[2],0,0,ierflg);
      min.mnparm(3,"HiggsB",pStart[3]+0.01,pStep[3],0,0,ierflg);
      
      
      min.mnexcm("MIGRAD",arglistS0,2,ierflg);
      min.mnexcm("MIGRAD",arglist,2,ierflg);
      if(ierflg!=0) {
	//std::cout << "ERROR: " << ierflg << std::endl;
	//std::cout << "S: " << params.S << std::endl;
	//min.SetPrintLevel(1);
	min.mnexcm("MIGRAD",arglist,2,ierflg);
	continue;
      }      
    }

    float thisL = getR(min);
    //std::cout <<thisL <<std::endl;
    if(TMath::IsNaN(thisL)) {
      std::cout << "NaN!  S=" << params.S << std::endl;
      continue;
    }

    //update minimum as needed
    if(thisL<curMin) {
      curMin=thisL;
      nGreater=0;
      bestS = params.S;
    } else {
      nGreater++;
    }

    if(nGreater>nReqGreater) break; //check nGreater
  }

  //transform to a -2 Delta LL
  
  //std::cout << curMin << std::endl;
  //std::cout << likeS0 << std::endl;
  
  //std::cout << sqrt(2*(TMath::Log(-1*curMin)-TMath::Log(-1*likeS0))) << "  " << bestS << std::endl;
  return std::make_pair(2*(TMath::Log(-1*curMin)-TMath::Log(-1*likeS0)) , bestS );

};


float getR(TMinuit& m) {
  //get results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  m.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  return amin;
}
