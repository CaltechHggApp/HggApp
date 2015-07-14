#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TLatex.h"

#include <iostream>

void drawAll_Better(TString folder,TString name,bool isData=true,TString smslbl="") {

  TFile file(folder+"/"+name+".root");

  //TList *l = file.GetListOfKeys();

  const int nCats=5;
  TString cats[nCats] = {"HighPt","Hbb","Zbb","HighRes","LowRes"};

  const int nHists=4;
  TString hists[nHists] = {"SignalRegion","SignalRegion_FineBin",
			   "SidebandRegion","SidebandRegion_FineBin"};

  TCanvas cv;
  cv.SetLogz();

  if(smslbl!="") {
    smslbl+="_";
  }

  for(int iCat=0; iCat<nCats; iCat++) {
    for(int iHist=0; iHist<nHists; iHist++) {
      if(!isData && iHist>1) break;
      TH2F* hist = (TH2F*)file.Get( Form("data_%s%s_%s",smslbl.Data(),cats[iCat].Data(),hists[iHist].Data()));
      hist->SetXTitle("M_{R} (GeV)");
      hist->SetYTitle("R^{2}");
      hist->Draw("COLZ");
      
      TLatex lbl0(0.1,0.96,TString("CMS Preliminary") + (isData ? "" : " Simulation"));
      lbl0.SetNDC();
      lbl0.SetTextSize(0.042);
      lbl0.Draw();
      
      TLatex lbl1(0.4,0.96,Form("%s Box",cats[iCat].Data()));
      lbl1.SetNDC();
      lbl1.SetTextSize(0.042);
      lbl1.Draw();
      
      TLatex lbl2(0.7,0.96,"#sqrt{s}= 8 TeV L = 19.78 fb^{-1}");
      lbl2.SetNDC();
      lbl2.SetTextSize(0.042);
      lbl2.Draw();
    
      cv.SaveAs(folder+"/figs/"+hist->GetName()+"_"+name+".png");    
      cv.SaveAs(folder+"/figs/"+hist->GetName()+"_"+name+".pdf");    
    }
  }

  

}
