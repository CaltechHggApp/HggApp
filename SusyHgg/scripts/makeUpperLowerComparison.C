#include "TFile.h"
#include "TH2F.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

#include "SigRegionBinning.h"

void makeUpperLowerComparison(TString dir,TString tag) {
  //data_LowRes_SidebandRegion_bkgShapeDown

  TFile f(dir+"/data.root");

  TCanvas cv;
  cv.SetLogy();

  for(int iBox=0;iBox<5;iBox++) {
    TString boxName = getRegionName( static_cast<SigRegionBinning::BinningRegion>(iBox) );
    TH2F* histDown = (TH2F*)f.Get( Form("data_%s_SidebandRegion_bkgShapeDown",boxName.Data()) );
    TH2F* histUp   = (TH2F*)f.Get( Form("data_%s_SidebandRegion_bkgShapeUp",boxName.Data()) );

    TH2F* histDown_Sig = SigRegionBinning::makeHistogram(static_cast<SigRegionBinning::BinningRegion>(iBox),"DownSideband");
    TH2F* histUp_Sig = SigRegionBinning::makeHistogram(static_cast<SigRegionBinning::BinningRegion>(iBox),"UpSideband");

    for(int iXbin=1; iXbin<=histDown->GetNbinsX(); iXbin++) {
      for(int iYbin=1; iYbin<=histDown->GetNbinsY(); iYbin++) {
	float x = histDown->GetXaxis()->GetBinCenter(iXbin);
	float y = histDown->GetYaxis()->GetBinCenter(iYbin);
	histDown_Sig->Fill(x,y, histDown->GetBinContent(iXbin,iYbin));
	histUp_Sig->Fill(x,y, histUp->GetBinContent(iXbin,iYbin));
      }
    }
    TH1D* histDown_1D = SigRegionBinning::make1DProj(histDown_Sig);
    TH1D* histUp_1D   = SigRegionBinning::make1DProj(histUp_Sig);

    histDown_1D->SetLineColor(kBlack);
    histUp_1D->SetLineColor(kRed);

    histDown_1D->SetLineWidth(2);
    histUp_1D->SetLineWidth(2);

    histDown_1D->SetMaximum(histDown_1D->GetMaximum()*10);

    histDown_1D->Draw();
    histUp_1D->Draw("SAME");
    
    TLegend leg(0.7,0.85,0.85,0.98);
    leg.SetFillColor(0);
    leg.SetBorderSize(0);

    leg.AddEntry(histDown_1D,"Lower Sideband", "l");
    leg.AddEntry(histUp_1D,"Upper Sideband", "l");
    leg.Draw("SAME");

    cv.SetLogy(1);
    cv.SaveAs(dir+Form("/figs/UpperLower_Comp_%s.png",boxName.Data()));
    cv.SaveAs(dir+Form("/figs/UpperLower_Comp_%s.pdf",boxName.Data()));

    TH1D* diff = (TH1D*)histDown_1D->Clone("diff");
    diff->Add(histUp_1D,-1);

    TH1D* error = (TH1D*)diff->Clone("ERROR");

    for(int i=1;i<=histDown_1D->GetNbinsX();i++) {
      float d = fabs(diff->GetBinContent(i));
      diff->SetBinContent(i, (d<0.001 ? 0.001:d) );
      error->SetBinError(i,TMath::Sqrt( histDown_1D->GetBinContent(i)+histUp_1D->GetBinContent(i))-0.01);
      error->SetBinContent(i,0.01);
    }

    error->SetMaximum(1E3);
    error->SetMinimum(0.01);
    error->SetFillColor(kRed);
    error->SetMarkerSize(0);
    diff->SetLineColor(kBlack);
    error->Draw("E2");
    diff->Draw("HSAME");
    TLegend leg2(0.55,0.85,0.95,0.99);
    leg2.SetFillColor(0);
    leg2.SetBorderSize(0);
    leg2.AddEntry(diff,"(Lower Sideband) - (Upper Sideband)","l");
    leg2.AddEntry(error,"1 #sigma Statistical Error","F");
    leg2.Draw("SAME");

    cv.SetLogy(1);
    cv.SaveAs(dir+Form("/figs/UpperLower_Diff_%s.png",boxName.Data()));
    cv.SaveAs(dir+Form("/figs/UpperLower_Diff_%s.pdf",boxName.Data()));

  }
}
