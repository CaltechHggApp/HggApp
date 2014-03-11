
#include "TH1D.h"

float getSigEff(const TH1D& hist, float mean) { //compute the sigma effective about the mean

  float integral = hist.Integral();

  int center_bin = hist.FindFixBin(mean);

  float N_bin=1;

  float target = 0.682689492;

  while(hist.Integral(center_bin-N_bin,center_bin+N_bin)/integral < target) N_bin++;

  return (N_bin-1)*hist.GetBinWidth(center_bin);
}


std::pair<float,float> getFastSIMScaleFactor(float eta, float pt) {
  float sfArray[6][5] = {
    {0.903,0.929,0.923,0.726,0.711},
    {0.949,0.949,0.954,0.803,0.810},
    {0.967,0.972,0.991,0.868,0.832},
    {0.974,0.979,1.021,0.896,0.862},
    {0.977,0.983,1.034,0.896,0.864},
    {0.972,0.977,0.998,0.893,0.854}
  };
  float sfError[6][5] = {
    {0.006,0.006,0.006,0.006,0.007},
    {0.002,0.002,0.002,0.002,0.003},
    {0.001,0.001,0.001,0.001,0.002},
    {0.001,0.001,0.001,0.001,0.001},
    {0.001,0.001,0.002,0.002,0.002},
    {0.002,0.002,0.002,0.002,0.003},    
  };

  int iEtaBin=0;
  int iPtBin =0;
  

  for(; iEtaBin < 5; iEtaBin++) {
    if(fabs(eta) <0.5*(iEtaBin+1)) break;
  }
  if(iEtaBin>4) iEtaBin=4;

  for(;iPtBin < 6; iPtBin++) {
    if( pt < (iPtBin+2)*10 ) break;
  }
  if(iPtBin >5) iPtBin=5;

  return std::make_pair( sfArray[iPtBin][iEtaBin], sfError[iPtBin][iEtaBin] );
}
