#include "TFile.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"

#include "SigRegionBinning.h"

#include <iostream>
#include <cstdio>

using namespace SigRegionBinning;

void MakeInvertedAnalysisTable(TString dir = "./") {
  TFile dataFile(dir+"/data.root");
  TFile SMHFile(dir+"/SMHiggs_SUM.root");
  TFile smsWH(dir+"/sms_ChiWH.root");
  TFile smsHH(dir+"/sms_ChiHH.root");

  for(int iBox=0;iBox<nBoxes;iBox++) {
    TString regName = getRegionName(static_cast<BinningRegion>(iBox));
    //data_75_200_LowRes_SignalRegion
    
    

    TH2F* sideband = (TH2F*)dataFile.Get(Form("data_%s_SidebandRegion",regName.Data()));
    TH2F* smhiggs = (TH2F*)SMHFile.Get(Form("data_%s_SignalRegion",regName.Data()));
    std::vector<TH2F*> smsWH_Hists;
    std::vector<TH2F*> smsHH_Hists;
    for(int m=125;m<=200;m+=25) {
      smsWH_Hists.push_back( (TH2F*)smsWH.Get(Form("data_0_%d_%s_SignalRegion",m,regName.Data())));
      smsHH_Hists.push_back( (TH2F*)smsHH.Get(Form("data_0_%d_%s_SignalRegion",m,regName.Data())));
    }

    std::cout << std::endl << regName << std::endl;
    const int RsqBin=2;
    std::cout << "R^2 cut:  " << sideband->GetYaxis()->GetBinLowEdge(RsqBin) << std::endl;
    std::cout << std::endl;
    for(int iBin=1;iBin<sideband->GetNbinsX();iBin++) {
      float bkg = sideband->Integral(iBin,-1,RsqBin,-1);
      float smh = smhiggs->Integral(iBin,-1,RsqBin,-1);
      if(bkg==0) continue;
      printf("% 3.0f & %0.1f & %0.1f & %0.1f ",sideband->GetXaxis()->GetBinLowEdge(iBin),bkg,smh,bkg+smh);
      for(int mm=0;mm<smsWH_Hists.size();mm++) {
	float sig = smsWH_Hists.at(mm)->Integral(iBin,-1,RsqBin,-1);
	printf(" & %0.1f & %0.3f",sig,sig/sqrt(bkg+smh));
      }
      std::cout << "\\\\" << std::endl;
    }

  }
}
