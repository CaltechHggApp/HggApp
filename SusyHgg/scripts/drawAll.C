#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TList.h"
#include "TRegexp.h"

#include <iostream>

void drawAll(TString folder,TString name,bool isData=true) {

  TFile file(folder+"/"+name+".root");

  TList *l = file.GetListOfKeys();

  TCanvas cv;
  cv.SetLogz();

  for(int i=0;i<l->GetEntries();i++) {
    TH1* hist = (TH1*)file.Get(l->At(i)->GetName());
    std::cout << hist->GetName() << std::endl;
    if( TString(hist->GetName()).Contains("Up") ||
	TString(hist->GetName()).Contains("Down") ) continue;

    if(TString(hist->GetName()).Contains("mgg_dist")){
      std::cout << "TH1D" << std::endl;
      hist->Rebin(50);

      hist->SetXTitle("m_{#gamma#gamma} [GeV]");
      hist->SetYTitle(Form("Events / %0.2f",hist->GetBinWidth(2)));
      hist->SetMarkerStyle(11);
      hist->SetMarkerSize(1.2);
      hist->Draw("P");
    }else{
      if(TString(hist->GetName()).Contains("SignalRegion") && isData) hist->SetMinimum(1);
      hist->SetXTitle("M_{R} [GeV]");
      hist->SetYTitle("R^{2}");
      hist->Draw("COLZ");
    }
    cv.SaveAs(folder+"/figs/"+hist->GetName()+"_"+name+".png");
  }

  

}
