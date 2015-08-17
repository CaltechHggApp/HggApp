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

#define nsamples 100000000
#define NBINS 42

void GetLEEpval( double localPval )
{
  double nsigma = TMath::Sqrt(2)*TMath::ErfcInverse( 2.*localPval );
  TRandom3 Gaus( 0 );
  int success = 0;
  for ( int i = 0; i < nsamples; i++ )
    {
      for( int j = 0; j < NBINS; j++ )
	{
	  if ( Gaus.Gaus() > nsigma )
	    {
	      success++;
	      break;
	    }
	}
    }
  std::cout << "success: " << success/(double)nsamples << std::endl;
  std::cout << "trial factor: " << (success/(double)nsamples)/localPval << std::endl;
  return;
}
