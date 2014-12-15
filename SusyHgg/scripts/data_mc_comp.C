#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"

TCanvas* buildDist(TFile* data_file, TFile* mc_file, TString reg,bool projX);

void data_mc_comp(TFile* data_file, TFile* mc_file,TString outputFolder) {
  const int nReg=5;
  TString regions[nReg] = {
    "data_HighPt_SidebandRegion",
    "data_Hbb_SidebandRegion",
    "data_Zbb_SidebandRegion",
    "data_HighRes_SidebandRegion",
    "data_LowRes_SidebandRegion",
  };

  for(int i=0;i<nReg;i++) {
    TCanvas *c = buildDist(data_file,mc_file,regions[i],true);
    c->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_MRdata.png");
    c->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_MRdata.pdf");
    c->SetLogy();
    c->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_MRdata_log.png");
    c->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_MRdata_log.pdf");
    delete c;
    TCanvas *d = buildDist(data_file,mc_file,regions[i],false);
    d->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_Rsqdata.png");
    d->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_Rsqdata.pdf");
    d->SetLogy();
    d->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_Rsqdata_log.png");
    d->SaveAs(outputFolder+"/dataMC_"+regions[i]+"_Rsqdata_log.pdf");

    delete d;
  }
  
}

TCanvas* buildDist(TFile* data_file, TFile* mc_file, TString reg,bool projX) {
  TCanvas *cv = new TCanvas();

  TH2F* data2d = (TH2F*)data_file->Get(reg);
  TH2F* mc2d = (TH2F*)mc_file->Get(reg);

  TH1D *data, *mc;

  if(projX) {
    data = data2d->ProjectionX("data");
    mc   = mc2d->ProjectionX("mc");
  }else{
    data = data2d->ProjectionY("data");
    mc   = mc2d->ProjectionY("mc");
  }

  mc->Scale(data->Integral()/mc->Integral());
  
  mc->SetFillColor(kRed);
  mc->SetLineColor(kRed);
  mc->SetFillStyle(3005);

  data->SetMarkerStyle(20);

  mc->SetMaximum(mc->GetMaximum()*1.3);
  mc->Draw();

  data->Draw("PSAME");

  TLegend leg(0.7,0.7,0.85,0.85);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);

  leg.AddEntry(data,"data","p");
  leg.AddEntry(mc,"MC","F");

  return cv;
}
