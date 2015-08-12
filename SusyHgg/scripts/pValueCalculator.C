#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TObject.h"
#include "TNamed.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TH1F.h"
#include "TRandom3.h"

#include "TBox.h"

#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooCurve.h"

#define nsamples 10000000
#define NBINS 100000

TH1F* GetPosteriorPdf( double pars[6] )
{
  double TF = pars[0];/*transfer factor*/
  double TF_Err = pars[1];/*trabsfer factor Unc.*/
  double Nside = pars[2];/*sideband statistics*/
  double Higgs = pars[3];/*higgs expected*/
  double HiggsErr = pars[4];/*higgs uncertainty*/
  double ShapeErr = pars[5];/*shape uncertainty*/

  /*ramdomize histo name*/
  TRandom3 histoName( 0 );
  //--------------------------------------
  //define a random var for each nuissance 
  //--------------------------------------
  TRandom3 gausTF( 0 );
  TRandom3 poissonNside( 0 );
  TRandom3 gausHiggs( 0 );
  TRandom3 gausShape( 0 );
  //---------------------------------
  //pdf for expected observation
  //---------------------------------
  TRandom3 poissonNobs( 0 );
  //histogram to store sampled Nobs
  int i_tmp = histoName.Integer( 100000 );
  TH1F* h_Nobs = new TH1F( Form("h_Nobs_%d", i_tmp), "h_Nobs", NBINS, -1, 4999. );
  for ( int i = 0; i < nsamples; i++ )
    {
      //--------------------------------
      //sample values for each parameter
      //--------------------------------
      double sTF    = gausTF.Gaus( TF, TF_Err );/*sampled value for TF*/
      double sNside = poissonNside.Poisson( Nside );
      double sHiggs = gausHiggs.Gaus( Higgs, HiggsErr );
      double sShape = gausShape.Gaus( 0.0, ShapeErr );
      
      //-----------------------
      //expected observation
      //using sampled values!!
      //-----------------------
      double poissonMean = sTF*sNside*( 1.0 + sShape ) + sHiggs;
      if (  poissonMean < .0 ) poissonMean = .0;
      double sNobs = poissonNobs.PoissonD( poissonMean );
      //double sNobs = poissonNobs.Gaus( poissonMean, TMath::Sqrt(poissonMean) );
      h_Nobs->Fill( sNobs );
    }
  
  return h_Nobs;
};

double getPval( double pars[6], int Nobs )
{
  TH1F* h_PDF = GetPosteriorPdf( pars );
  double p_val = -99;
  int bin = h_PDF->FindBin( Nobs );
  if ( h_PDF->GetMean() < Nobs )
    {
      //double delta = Nobs - h_PDF->GetMean();
      //int lbin = h_PDF->FindBin( h_PDF->GetMean() - delta );
      //p_val =  0.5*( h_PDF->Integral( bin, NBINS ) +  h_PDF->Integral( 1, lbin ) )/h_PDF->Integral();
      p_val = h_PDF->Integral( bin, NBINS )/h_PDF->Integral();
    }
  else
    {
      //double delta = h_PDF->GetMean() - Nobs;
      //int hbin = h_PDF->FindBin( h_PDF->GetMean() + delta );
      //p_val =  0.5*( h_PDF->Integral( 1, bin ) + h_PDF->Integral( hbin, NBINS ) )/h_PDF->Integral();
      p_val = h_PDF->Integral( 1, bin )/h_PDF->Integral();
    }
  
  return p_val;
};

double getNsigma( double pars[6], int Nobs )
{
  return TMath::Sqrt(2)*TMath::ErfcInverse( 2.*getPval( pars, Nobs ) );
};

double getNsigma( double pval )
{
  return TMath::Sqrt(2)*TMath::ErfcInverse( 2.*pval );
};
