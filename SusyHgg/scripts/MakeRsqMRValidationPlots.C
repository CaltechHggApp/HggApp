
#include "TH2F.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TFile.h"

#include "SelectionStrings.h"
#include "SigRegionBinning.h"

#include <iostream>

enum sel_t {kReal,kShift,kSingleIso};

void makeValidationPlots(TTree* tree, sel_t type,bool isMC) {

  SelectionStrings strings;

  TString tag="";

  TString lowString="",midString="",highString="";

  switch(type) {
  case kReal:
    tag = "_realselection";
    break;
  case kShift:
    tag = "_side";
    break;
  case kSingleIso:
    tag = "_singleIso";
    break;
  }
  if(isMC) {
    tag+="__DiPhotonJets";
  }
  Double_t Red[] = {0.00, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00};
  Double_t Green[] ={0.00, 0.70, 0.90, 1.00, 0.90, 0.70, 0.00};
  Double_t Blue[] = {1.00, 1.00, 1.00, 1.00, 0.90, 0.70, 0.00};
  Double_t Length[] =  {0.00, 0.20, 0.35, 0.50, 0.65, 0.8, 1.00};

  TCanvas cv;

  std::vector<TH2F*> nSigs;
  std::vector<TH2F*> nSigs_gauss;

  cv.SetLogx();
  for(int iBox=0; iBox<SigRegionBinning::nBoxes; iBox++) {
    TString boxName = SigRegionBinning::getRegionName(static_cast<SigRegionBinning::BinningRegion>(iBox));
    switch(type) {
    case kReal:
      lowString = strings.baseSelection+" && mgg>103 && mgg<120";
      midString = strings.baseSelection+" && "+strings.mggSigRegion[iBox];
      highString = strings.baseSelection+" && mgg>131 && mgg<160";
      break;
    case kShift:
      lowString = strings.baseSelection+" && mgg>130 && mgg<140";
      midString = strings.baseSelection+" && mgg>140 && mgg<150";
      highString = strings.baseSelection+" && mgg>150 && mgg<160";
      break;
    case kSingleIso:
      lowString = TString("(pho1_pt>40 && pho2_pt>25 && abs(pho1_eta)<1.48 && abs(pho2_eta)<1.48 && (pho1_pass_iso || pho2_pass_iso) && !(pho1_pass_iso && pho2_pass_iso))")+" && mgg>103 && mgg<120";
      midString = TString("(pho1_pt>40 && pho2_pt>25 && abs(pho1_eta)<1.48 && abs(pho2_eta)<1.48 && (pho1_pass_iso || pho2_pass_iso) && !(pho1_pass_iso && pho2_pass_iso))")+" && "+strings.mggSigRegion[iBox];
      highString = TString("(pho1_pt>40 && pho2_pt>25 && abs(pho1_eta)<1.48 && abs(pho2_eta)<1.48 && (pho1_pass_iso || pho2_pass_iso) && !(pho1_pass_iso && pho2_pass_iso))")+" && mgg>131 && mgg<160";      
      break;
    }


    TH2F* low = SigRegionBinning::makeHistogram(static_cast<SigRegionBinning::BinningRegion>(iBox),"low");
    TH2F* mid = SigRegionBinning::makeHistogram(static_cast<SigRegionBinning::BinningRegion>(iBox),"mid");
    TH2F* high = SigRegionBinning::makeHistogram(static_cast<SigRegionBinning::BinningRegion>(iBox),"high");
  
    TH2F* nsig = SigRegionBinning::makeHistogram(static_cast<SigRegionBinning::BinningRegion>(iBox),"nsig_"+boxName);
    TH2F* nsig_gauss = SigRegionBinning::makeHistogram(static_cast<SigRegionBinning::BinningRegion>(iBox),"nsig_gauss_"+boxName);

    low->SetMinimum(0.1);
    mid->SetMinimum(0.1);
    high->SetMinimum(0.1);

    tree->Project("low","Rsq:MR",lowString+" && "+strings.boxDefs[iBox]);
    tree->Project("mid","Rsq:MR",midString+" && "+strings.boxDefs[iBox]);
    tree->Project("high","Rsq:MR",highString+" && "+strings.boxDefs[iBox]);

    //SigRegionBinning::formatSigRegionPlot(low);
    //SigRegionBinning::formatSigRegionPlot(mid);
    //SigRegionBinning::formatSigRegionPlot(high);

    //get the difference in the prediction
    TH2F* low_norm = (TH2F*)low->Clone("low_norm");
    TH2F* high_norm = (TH2F*)high->Clone("high_norm");

    low_norm->Scale(1./low_norm->Integral());
    high_norm->Scale(1./high_norm->Integral());

    TH2F* pred_diff = (TH2F*)high_norm->Clone("pred_diff");
    pred_diff->Add(low_norm,-1);

    TH2F* norm_av = (TH2F*)low_norm->Clone("norm_av");
    norm_av->Add(high_norm);
    norm_av->Scale(0.5);

    TH2F* pred_per_diff = (TH2F*)pred_diff->Clone("pred_per_diff");
    pred_per_diff->Divide(norm_av);
    cv.SetLogx(1);
    cv.SetLogy(0);
    pred_per_diff->Draw("COLZ TEXT");
    cv.SaveAs("RsqMR_high_minus_low_div_av_"+boxName+tag+".png");
    cv.SetLogx(1);
    cv.SetLogy(0);

    TH1D* low_1D = SigRegionBinning::make1DProj(low);
    TH1D* high_1D = SigRegionBinning::make1DProj(high);
    
    pred_per_diff->SetYTitle("high - low / (high+low/2)");
    TH1D* pred_per_diff_1D = SigRegionBinning::make1DProj(pred_per_diff);
    cv.SetLogx(0);
    cv.SetLogy(0);
    pred_per_diff_1D->Draw();
    cv.SaveAs("RsqMR_1D_high_minus_low_div_av_"+boxName+tag+".png");
    cv.SetLogx(1);
    cv.SetLogy(0);
    
    TFile out("RsqMR_1D_high_minus_low_div_av_"+boxName+tag+".root","RECREATE");
    low->Write();
    high->Write();
    low_1D->Write();
    high_1D->Write();
    pred_per_diff_1D->Write();
    pred_per_diff->Write();
    out.Close();



    TH2F* sideband_tot = (TH2F*)low->Clone("sideband_tot");
    sideband_tot->Add(high);


    sideband_tot->Draw("COLZ TEXT");
    cv.SetLogz();
    cv.SaveAs("RsqMR_low_plus_high_"+boxName+tag+".png");
    cv.SetLogz(0);

    float sf = mid->Integral()/sideband_tot->Integral();
    std::cout << "sf: " << sf << std::endl;
    TH2F* sideband_pred = (TH2F*)sideband_tot->Clone("sideband_pred");

    sideband_pred->Scale(sf);

    sideband_pred->Draw("COLZ TEXT");
    cv.SetLogz();
    cv.SaveAs("RsqMR_sideband_pred_"+boxName+tag+".png");
    cv.SetLogz(0);

    mid->Draw("COLZ TEXT");
    cv.SetLogz();
    cv.SaveAs("RsqMR_mid_"+boxName+tag+".png");
    cv.SetLogz(0);

    TH1D* pred_1D = SigRegionBinning::make1DProj(sideband_pred,true,sf);
    TH1D* mid_1D  = SigRegionBinning::make1DProj(mid);
    cv.SetLogx(0);
    cv.SetLogy();
    mid_1D->SetMarkerStyle(20);
    mid_1D->SetMarkerColor(kBlack);
    mid_1D->SetMarkerSize(1.4);

    pred_1D->SetFillColor(kRed);
    pred_1D->SetFillStyle(3001);
    pred_1D->SetLineColor(kRed);
    pred_1D->SetLineWidth(2);
    pred_1D->SetMarkerSize(0);
    pred_1D->Draw("L E2");
    pred_1D->Draw("LSAME");
    mid_1D->Draw("PSAME");
    cv.SaveAs("RsqMR_1D_obs_sideband_pred_"+boxName+tag+".png");    
    cv.SetLogx(1);
    cv.SetLogy(0);


    TH2F* perdiff = (TH2F*)mid->Clone("perdiff");
    perdiff->Add(sideband_pred,-1);
    perdiff->Divide(perdiff);
    perdiff->Draw("COLZ TEXT");

    cv.SaveAs("RsqMR_percentDiff_"+boxName+tag+".png");
    

    for(int iXbin=1; iXbin<=sideband_pred->GetNbinsX(); iXbin++) {
      for(int iYbin=1; iYbin<=sideband_pred->GetNbinsY(); iYbin++) {
	float obs = mid->GetBinContent(iXbin,iYbin);
	float exp = sideband_pred->GetBinContent(iXbin,iYbin);
	float err = TMath::Sqrt(sideband_tot->GetBinContent(iXbin,iYbin))*sf;
	float ns = fabs(TMath::NormQuantile(SigRegionBinning::pValue(obs,exp,err/exp)/2));
	if(obs<exp) ns*=-1;

	nsig->SetBinContent(iXbin,iYbin,ns);
	std::cout << "\t" << iXbin << "  " << iYbin << "  " << nsig->GetBinContent(iXbin,iYbin) << std::endl;

	float gauss_err = TMath::Sqrt(TMath::Power(err,2)+exp);
	if(gauss_err==0) gauss_err=1;
	nsig_gauss->SetBinContent(iXbin,iYbin,(obs-exp)/gauss_err);

      }
    }
    
    SigRegionBinning::formatSigRegionPlot(nsig); 
    SigRegionBinning::formatSigRegionPlot(nsig_gauss); 
    //cv.SetLogx();

    nSigs.push_back(nsig);
    nSigs_gauss.push_back(nsig_gauss);
    

    delete low;
    delete mid;
    delete high;
    delete sideband_tot;
    delete sideband_pred;

  }

  for(int iBox=0; iBox<SigRegionBinning::nBoxes; iBox++) {
    TString boxName = SigRegionBinning::getRegionName(static_cast<SigRegionBinning::BinningRegion>(iBox));
    TColor::CreateGradientColorTable(7,Length,Red,Green,Blue,999);
    TH2F* nsig = nSigs.at(iBox);

    nsig->SetMaximum(5.1);
    nsig->SetMinimum(-5.1);
    nsig->SetContour(999);


    nsig->Draw("COLZ TEXT0");
    cv.SaveAs("RsqMR_low_plus_high_RAWSUM_NSIGMA_"+boxName+tag+".png");
    
    TH2F* nsig_gauss = nSigs_gauss.at(iBox);
    nsig_gauss->SetMaximum(5.1);
    nsig_gauss->SetMinimum(-5.1);
    nsig_gauss->SetContour(999);

    nsig_gauss->Draw("COLZ");
    cv.SaveAs("RsqMR_low_plus_high_RAWSUM_NSIGMAGAUSS_"+boxName+tag+".png");
    

    TH1D* nsig_gauss_1D = SigRegionBinning::make1DProj(nsig_gauss);

    nsig_gauss_1D->SetLineColor(kBlack);
    nsig_gauss_1D->SetLineWidth(1.5);
    nsig_gauss_1D->SetMinimum(-5.1);
    nsig_gauss_1D->SetMaximum(5.1);
    nsig_gauss_1D->SetYTitle("Number of Sigma");
    cv.SetLogx(0);
    cv.SetLogy(0);
    nsig_gauss_1D->Draw();
    cv.SaveAs("RsqMR_low_plus_high_1D_RAWSUM_NSIGMAGAUSS_"+boxName+tag+".png");
    cv.SetLogx(1);
  }
}
