#include "TMinuit.h"
#include "TMath.h"

#include <iostream>

#include "TH1D.h"

//the parameters of the bin we are in
struct par {
  double obs=2;
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
  f = 1;
  f *= TMath::Gaus(par[0],params.UpperSF,params.UpperSFe); // upper SF compatibility
  f *= TMath::Gaus(par[1],params.LowerSF,params.LowerSFe); // lower SF compatibility
  f *= TMath::PoissonI(params.Nupper,par[2]);                    // upper sideband compatibility
  f *= TMath::PoissonI(params.Nlower,par[2]*par[0]/par[1]);      // lower sideband compatibility
  f *= TMath::Gaus(par[3],params.Higgs,params.HiggsErr);        // Higgs background compatibility
  f *= TMath::PoissonI(params.obs,params.S+par[0]*par[2]+par[3]);            // S compatibility
  f*=-1;
}


TH1D* profileMinuit(float max) {
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



  TH1D* output = new TH1D("output","",max/100,0,max);
  for(int i=0;i<output->GetNbinsX(); i++) {
    min.mnparm(0,"UpperSF",pStart[0],pStep[0],0,0,ierflg);
    min.mnparm(1,"LowerSF",pStart[1],pStep[1],0,0,ierflg);
    min.mnparm(2,"UpperB",pStart[2],pStep[2],0,0,ierflg);
    min.mnparm(3,"HiggsB",pStart[3],pStep[3],0,0,ierflg);

    params.S = output->GetBinLowEdge(i+1);
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
	 std::cout << "ERROR: " << ierflg << std::endl;
	 std::cout << "S: " << params.S << std::endl;
	 //min.SetPrintLevel(1);
	 min.mnexcm("MIGRAD",arglist,2,ierflg);
	continue;
       }
    }
    
    //get results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    min.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    std::cout.precision(10);
    //std::cout <<  amin << std::endl;
    output->SetBinContent(i+1,-1*amin);
  }
//   std::cout << "amin: " << amin << std::endl;
//   std::cout << "edm: "  << edm  << std::endl;
//   std::cout << "errdef: "  << errdef  << std::endl;
//   std::cout << "nvpar: "  << nvpar  << std::endl;
//   std::cout << "nparx: "  << nparx  << std::endl;
//   std::cout << "icstat: "  << icstat  << std::endl;

  //transform to a -2 Delta LL


  double maximum = output->GetMaximum();
  params.S = 9.3;
  //params.S = output->GetBinLowEdge(output->GetMaximumBin());
  //std::cout << "S: " << params.S << std::endl;
  //min.SetPrintLevel(1);
  //min.mnexcm("MIGRAD",arglist,2,ierflg);
  
  for(int i=1;i<=output->GetNbinsX();i++) {
    //std::cout << 2*(TMath::Log(maximum)-TMath::Log(output->GetBinContent(i))) << std::endl;
    output->SetBinContent(i,2*(TMath::Log(maximum)-TMath::Log(output->GetBinContent(i))));
  }
  
  return output;
};
