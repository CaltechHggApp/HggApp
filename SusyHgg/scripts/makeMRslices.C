#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"

#include "TF1.h"

void getSlice(TFile *f,TString reg,TString name,int col, std::vector<TH1D*>* hists) {
  TH2F* h2d = (TH2F*)f->Get(reg);
  TH1D* r0p00 = h2d->ProjectionX(name+"_r0p00");
  TH1D* r0p05 = h2d->ProjectionX(name+"_r0p05",2);
  TH1D* r0p10 = h2d->ProjectionX(name+"_r0p10",3);

  r0p00->SetFillColor(col);
  r0p05->SetFillColor(col);
  r0p10->SetFillColor(col);

  r0p05->SetAxisRange(1E-1,800,"Y");
  r0p05->SetAxisRange(0,1200,"X");

  r0p10->SetAxisRange(1E-1,800,"Y");
  r0p10->SetAxisRange(0,1200,"X");

  hists->push_back(r0p00);
  hists->push_back(r0p05);
  hists->push_back(r0p10);

}

void makeMRslices(TFile *data, TFile *fggH,TFile *fvbfH, TFile *fwzH, TFile *fttH,TString outputDir) {
  std::vector<TH1D*> obs;
  std::vector<TH1D*> bkg;
  std::vector<TH1D*> ggH;
  std::vector<TH1D*> vbfH;
  std::vector<TH1D*> wzH;
  std::vector<TH1D*> ttH;

  std::vector<float> rVals = {0.00,0.05,0.10};

  getSlice(data,"data_HighRes_SignalRegion","obs",1,&obs);
  getSlice(data,"data_HighRes_SidebandRegion","bkg",2,&bkg);
  getSlice(fggH,"data_HighRes_SignalRegion","ggH",3,&ggH);
  getSlice(fvbfH,"data_HighRes_SignalRegion","vbfH",4,&vbfH);
  getSlice(fwzH,"data_HighRes_SignalRegion","wzH",5,&wzH);
  getSlice(fttH,"data_HighRes_SignalRegion","ttH",6,&ttH);
  

  for(int i=0;i<obs.size();i++) {

    obs.at(i)->SetMarkerStyle(20);

    THStack rStack("rStack","");

    rStack.Add(ggH.at(i));
    rStack.Add(vbfH.at(i));
    rStack.Add(wzH.at(i));
    rStack.Add(ttH.at(i));
    rStack.Add(bkg.at(i));
    
    TLegend leg(0.5,0.6,0.85,0.9);
    leg.SetFillColor(0);
    leg.SetBorderSize(0);
    
    leg.AddEntry(obs.at(i),"signal region data","p");
    leg.AddEntry(bkg.at(i),"sideband background prediction","F");
    leg.AddEntry(ggH.at(i),"pp#rightarrow H","F");
    leg.AddEntry(vbfH.at(i),"pp#rightarrow qqH","F");
    leg.AddEntry(wzH.at(i),"pp#rightarrow VH","F");
    leg.AddEntry(ttH.at(i),"pp#rightarrow ttH","F");

    TCanvas cv;
    cv.SetLogy();

    obs.at(i)->Draw("P");
    rStack.Draw("SAME");
    obs.at(i)->Draw("PSAME");
    leg.SetHeader(Form("R^{2} > %0.2f",rVals.at(i)));
    leg.Draw();

    cv.SaveAs(outputDir+Form("/mrDist_obs_bkg_Rsq_0p%02d.png",int(100*rVals.at(i))));
    cv.SaveAs(outputDir+Form("/mrDist_obs_bkg_Rsq_0p%02d.pdf",int(100*rVals.at(i))));

    TH1D* total = (TH1D*) bkg.at(i)->Clone("total_bkg");
    std::cout << total->Integral() << std::endl;
    total->Add(ggH.at(i));
    total->Add(vbfH.at(i));
    total->Add(wzH.at(i));
    total->Add(ttH.at(i));
    std::cout << total->Integral() << std::endl;


    TF1 bkgFit("bkgFit","expo(0)",150,500);
    //bkgFit.SetParameter(0,TMath::Log(total->GetMaximum()));
    //bkgFit.SetParameter(1,-1);
    
    //bkgFit.SetParLimits(0,1,1E9);
    //bkgFit.SetParLimits(1,-100,-0.00001);

    obs.at(i)->Fit(&bkgFit,"MR");

    for(int j=1;j<total->GetNbinsX()+1;j++) {
      total->SetBinContent(j, total->GetBinContent(j)*50/total->GetBinWidth(j));
      obs.at(i)->SetBinContent(j, obs.at(i)->GetBinContent(j)*50/total->GetBinWidth(j));
    }


    TCanvas cv2;
    cv2.SetLogy();
    total->SetFillColor(kRed);
    
    total->Draw("H");
    bkgFit.SetLineColor(kGreen);
    bkgFit.Draw("LSAME");
    obs.at(i)->Draw("PSAME");
    TFile out(outputDir+Form("/mrDist_obs_bkg_fit_Rsq_0p%02d.root",int(100*rVals.at(i))),"RECREATE");
    total->Write();
    bkgFit.Write();
    obs.at(i)->Write();
    out.Close();

    cv2.SaveAs(outputDir+Form("/mrDist_obs_bkg_fit_Rsq_0p%02d.png",int(100*rVals.at(i))));
 }
}

void makeMRslices_thisDir() {
  makeMRslices(TFile::Open("data.root"),TFile::Open("ggH.root"),TFile::Open("vbfH.root"),TFile::Open("wzH.root"),TFile::Open("ttH.root"),"figs/");
}


