#include "TMinuit.h"
#include "TMath.h"

#include <iostream>

#include "TH1D.h"
#include "TRandom.h"

//the parameters of the bin we are in
struct par {
  double obs=5;
  double Nside=3;
  double sf=0.162,sfe=0.0016;
  double combAddErr=0.5; //additional % error on the combinatorial bkg
  double Higgs=0.05,HiggsErr=0.012;
  double S=0;
} params;


double gaus(double x, double m, double s) {
  return 1/sqrt(TMath::Pi()*s)*TMath::Exp( -1*(x-m)*(x-m)/s/s);
}

double pois(double N, double alpha) {
  return TMath::Power(alpha,N)*TMath::Exp(-1*alpha);
}


//nuisance parameters
//float SF,B,Higgs;

// npar should be 4 (free) parameters
// SF,B,HiggsB,bkgMMult
void negLikelihood(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t flag) {
  if(params.S+par[0]*par[2]+par[3] <0) { 
    f=0; //truncate so we can scan negative S
    return;
  }
  f = 1;
  f *= TMath::Gaus(par[0],params.sf,params.sfe);//SF
  //f *= TMath::Poisson(params.Nside,par[1]);//sideband compatibility
  f *= TMath::PoissonI(params.Nside,par[1]);
  f *= TMath::Gaus(par[2],params.Higgs,params.HiggsErr);//Higgs background compatibility
  f *= TMath::Gaus(par[3],0,params.combAddErr);//background shape compatibility
  //f *= TMath::Poisson(params.obs,params.S+par[0]*par[1]*(1+par[3])+par[2]);//signal+bkg compatibility
  f *= TMath::PoissonI(params.obs,params.S+par[0]*par[1]*(1+par[3])+par[2]);
  f *= -1;
}

//Allow last poison to be float instead of integer
void negLikelihoodMod(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t flag) {
  if(params.S+par[0]*par[2]+par[3] <0) {
    f=0; //truncate so we can scan negative S                                                                           
    return;
  }
  f = 1;
  f *= TMath::Gaus(par[0],params.sf,params.sfe);//SF                                                    
  f *= TMath::Poisson(params.Nside,par[1]);//sideband compatibility
  f *= TMath::Gaus(par[2],params.Higgs,params.HiggsErr);//Higgs background compatibility       
  f *= TMath::Gaus(par[3],0,params.combAddErr);//background shape compatibility       
  f *= TMath::Poisson(params.obs,par[0]*par[1]*(1+par[3])+par[2]);//no signal, bkg compatibiblity
  f *= -1.0;
}

float getR(TMinuit& m);

//just return the -2delta LL , best S
std::pair<float,float> profileSimpleNoHistMinuit( float step = 0.01, bool _singlePoint = false, double signal = .0,
						  bool _profileNobs = false, double obs = 0.0 ) {
  //    std::cout << std::endl;
  //    std::cout << params.obs << " " << params.Nupper << " " << params.Nlower << " " << params.UpperSF << " " << params.UpperSFe << " " << std::endl;
  //    std::cout << params.LowerSF << " " << params.LowerSFe << " " << params.Higgs << " " << params.HiggsErr << " " <<  std::endl;
  
  const int nPar = 4;
  
  TMinuit min(nPar);
  Int_t ierflg = 0;
  
  if ( !_profileNobs )
    {
      min.SetFCN( negLikelihood );
    }
  else
    {
      min.SetFCN( negLikelihoodMod );
    }
  
  
  double pStart[nPar]={params.sf,params.Nside,params.Higgs,0};
  double pStep [nPar]={0.0001,0.0001,0.0001,0.0001};
  
  
  Double_t arglist[10];
  Double_t arglistS0[10];
  
  arglist[0] = 5000;
  arglist[1] = 2.;
  
  arglistS0[0] = 5000;
  arglistS0[1] = 0.;
  
  min.SetPrintLevel(-1);
  
  float curMin = 99999999999; // store the minimum value of the likelihood found
  float bestS = 0;
  const int nReqGreater = int(1./step); // number of point greater than the curMin required to declare global min
  int nGreater = 0;          // number of points currently found greater than curMin
  
  float likeS0 = 0; //store the likelihood at S=0
  
  //do the S=0 likelihood
  
  const int nMax = 1E6; // never exceed this many calls
  
  params.S = 0;
  min.mnparm(0,"SF",pStart[0],pStep[0],0,0,ierflg);
  min.mnparm(1,"B",pStart[1],pStep[1],0,0,ierflg);
  min.mnparm(2,"HiggsB",pStart[2],pStep[2],0,0,ierflg);
  min.mnparm(3,"BkgShape",pStart[3],pStep[3],0,0,ierflg);
  
  //min.SetPrintLevel(1);
  
  min.mnexcm("MIGRAD",arglistS0,2,ierflg);
  min.mnexcm("MIGRAD",arglist,2,ierflg);
  
  likeS0 = getR(min);
  
  if ( _singlePoint )
    {
      params.S = signal;
      min.mnparm(0,"SF",pStart[0],pStep[0],0,0,ierflg);
      min.mnparm(1,"B",pStart[1],pStep[1],0,0,ierflg);
      min.mnparm(2,"HiggsB",pStart[2],pStep[2],0,0,ierflg);
      min.mnparm(3,"BkgShape",pStart[3],pStep[3],0,0,ierflg);
      
      if( _profileNobs ) params.obs = obs;
      //min.SetPrintLevel(1);
      min.mnexcm("MIGRAD",arglistS0,2,ierflg);
      min.mnexcm("MIGRAD",arglist,2,ierflg);
      
      
      if( TMath::IsNaN( getR(min) ) ) 
	{
	  std::cout << "NaN!  obs=" << params.obs << std::endl;
	  return std::make_pair( 0, 0 );
	}
      return std::make_pair( 2*TMath::Log( -1*getR(min) ), signal );
      //return std::make_pair( -1*getR(min), signal );
    }
  
  if(TMath::IsNaN(likeS0)) {
    std::cout << "FATAL ERROR: S=0 likelihood is NaN" << std::endl;
    return std::make_pair(0.,0.);
  }
  min.SetPrintLevel(-1);
  
  curMin = likeS0;
  
  //scan S
  int nRun=0;
  
  while(true) {
    if(nRun++ > nMax) return std::make_pair(0,0);
    params.S+=step;
    min.mnparm(0,"SF",pStart[0],pStep[0],0,0,ierflg);
    min.mnparm(1,"B",pStart[1],pStep[1],0,0,ierflg);
    min.mnparm(2,"HiggsB",pStart[2],pStep[2],0,0,ierflg);
    min.mnparm(3,"BkgShape",pStart[3],pStep[3],0,0,ierflg);
    
    
    //std::cout << params.S << std::endl;
    min.mnexcm("MIGRAD",arglistS0,2,ierflg);
    min.mnexcm("MIGRAD",arglist,2,ierflg);
    if(ierflg!=0) {
      min.mnparm(0,"SF",pStart[0],pStep[0],0,0,ierflg);
      min.mnparm(1,"B",pStart[1],pStep[1],0,0,ierflg);
      min.mnparm(2,"HiggsB",pStart[2],pStep[2],0,0,ierflg);
      min.mnparm(3,"BkgShape",pStart[3],pStep[3],0,0,ierflg);
      
      
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
    min.mnparm(0,"SF",pStart[0],pStep[0],0,0,ierflg);
    min.mnparm(1,"B",pStart[1],pStep[1],0,0,ierflg);
    min.mnparm(2,"HiggsB",pStart[2],pStep[2],0,0,ierflg);
    min.mnparm(3,"BkgShape",pStart[3],pStep[3],0,0,ierflg);
    
    
    //std::cout << params.S << std::endl;
    min.mnexcm("MIGRAD",arglistS0,2,ierflg);
    min.mnexcm("MIGRAD",arglist,2,ierflg);
    if(ierflg!=0) {
      min.mnparm(0,"SF",pStart[0],pStep[0],0,0,ierflg);
      min.mnparm(1,"B",pStart[1],pStep[1],0,0,ierflg);
      min.mnparm(2,"HiggsB",pStart[2],pStep[2],0,0,ierflg);
      min.mnparm(3,"BkgShape",pStart[3],pStep[3],0,0,ierflg);
      
      
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
};

TH1D* getTwoLogLikelihood( bool _profileNobs = false )
{
  float step = 0.01;
  double x_h = 40., x_l = 0.0;
  if ( _profileNobs ) x_h = 3.0*params.Nside*params.sf;
  int nbins = (int) ( x_h - x_l )/step;
  double _nll[nbins];
  double _min_nll = 999999999.;
  TRandom a(0);
  TString _hname = Form( "_h_%d", a.Integer(10000) );
  TH1D* _h = new TH1D( _hname, "LogLikelihood", nbins, x_l, x_h );
  for ( int i = 1; i <= nbins; i++ )
    {
      double signal = x_l + step*( i - 1 );
      double ll;
      if ( !_profileNobs )
	{
	  ll = -1.0*profileSimpleNoHistMinuit( step, true, signal ).first;
	}
      else
	{
	  ll = -1.0*profileSimpleNoHistMinuit( step, true, 0.0, true, signal ).first;
	}
      if ( ll < _min_nll ) _min_nll = ll;//store nll minimum
      _nll[i-1] = ll;
      //_h->SetBinContent( i, ll );
    }
  for ( int i = 1; i <= nbins; i++ )
    {
      _h->SetBinContent( i, _nll[i-1] - _min_nll ); 
    }
  return _h;
};

