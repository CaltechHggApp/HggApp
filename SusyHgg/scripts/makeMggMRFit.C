#include "TTree.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooExtendPdf.h"
#include "RooStats/SPlot.h"
#include "RooStats/ModelConfig.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"

#include "TLatex.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

#include <iostream>

TString makeDoubleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeSingleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeSinglePow(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeDoublePow(TString tag, RooRealVar& mgg,RooWorkspace& w);


RooWorkspace* makeMggMRFit(TTree* tree,bool doBkg=true, float forceSigma=-1,float forceMu=-1) {
  RooWorkspace *ws = new RooWorkspace("ws","");
  RooRealVar mgg("mgg","m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(40);

  mgg.setRange("sideband_low", 103,120);
  mgg.setRange("sideband_high",131,160);

  RooRealVar MR("MR","",200,3000,"GeV");
  MR.setBins(60);
  RooRealVar Rsq("Rsq","",0,1,"GeV");
  Rsq.setBins(20);
  RooRealVar hem1_M("hem1_M","",0,2000,"GeV");
  hem1_M.setBins(40);
  RooRealVar hem2_M("hem2_M","",0,2000,"GeV");
  hem2_M.setBins(40);
  RooRealVar ptgg("ptgg","p_{T}^{#gamma#gamma}",0,500,"GeV");
  ptgg.setBins(50);

  RooDataSet data("data","",tree,RooArgSet(mgg,MR,Rsq,hem1_M,hem2_M,ptgg));

  TString tag_mgg = makeDoubleExp("mgg_bkg",mgg,*ws);
  TString tag_MR = makeModExp("MR_bkg",MR,*ws);

  RooAbsPdf * mggBkg = ws->pdf("mgg_bkg_"+tag_mgg);
  RooAbsPdf * mrBkg = ws->pdf("MR_bkg_"+tag_MR);

  if(mggBkg==0 || mrBkg==0) {
    std::cout << mggBkg << "\n" << mrBkg << std::endl;
    return 0;
  }
  
  RooProdPdf bkgModel("bkgModel","",*mggBkg,*mrBkg);

  RooRealVar mu_MR("mu_MR","",500,200,2000);
  mu_MR.setVal(500.);
  mu_MR.setConstant(true);
  RooRealVar mu_mgg("mu_mgg","",125,120,130);
  if(forceMu>0) {
    mu_mgg.setRange(110,170);
    mu_mgg.setVal(forceMu);
    mu_mgg.setConstant(true);
  }
  RooRealVar sigma_mgg("sigma_mgg","",1,0.2,6);
  if(forceSigma>0) {
    sigma_mgg.setVal(forceSigma);
    sigma_mgg.setConstant(true);
  }

  RooRealVar fsigma_MR("fsigma_MR","",0.3,0.1,0.6);
  RooFormulaVar sigma_MR("sigma_MR","","@0*@1",RooArgList(mu_MR,fsigma_MR));

  RooGaussian sigModel_mgg("sigModel_mgg","",mgg,mu_mgg,sigma_mgg);
  RooGaussian sigModel_MR("sigModel_MR","",MR,mu_MR,sigma_MR);

  RooProdPdf sigModel("sigModel","",sigModel_mgg,sigModel_MR);

  RooRealVar Nbkg("Nbkg","",10,1,1E9);
  RooRealVar Nsig("Nsig","",1,0,1E4);  


  RooAddPdf fitModel("fitModel","",RooArgList(bkgModel,sigModel),RooArgList(Nbkg,Nsig));

  fitModel.fitTo(data,RooFit::Strategy(0));
  RooFitResult*res = fitModel.fitTo(data,RooFit::Strategy(2),RooFit::Save(true),RooFit::Extended(kTRUE));
  res->SetName("fitres");

  RooPlot *fmgg = mgg.frame();
  fmgg->SetName("mgg_frame");
  data.plotOn(fmgg);
  fitModel.plotOn(fmgg);

  TLatex muL(0.62,0.8,Form("#mu = %0.1f #pm %0.2f GeV",mu_mgg.getVal(),mu_mgg.getError()));
  muL.SetTextSize(0.042);
  muL.SetNDC();

  TLatex sigL(0.62,0.72,Form("#sigma = %0.1f #pm %0.2f GeV",sigma_mgg.getVal(),sigma_mgg.getError()));
  sigL.SetTextSize(0.042);
  sigL.SetNDC();

  TLatex nL(0.62,0.64,Form("N_{sig} = %0.1f #pm %0.2f",Nsig.getVal(),Nsig.getError()));
  nL.SetTextSize(0.042);
  nL.SetNDC();

  fmgg->addObject(&muL);
  fmgg->addObject(&sigL);
  fmgg->addObject(&nL);
  
  RooPlot *fMR = MR.frame();
  fMR->SetName("MR_frame");
  data.plotOn(fMR);
  fitModel.plotOn(fMR);

  /*
  RooStats::SPlot* sdata = new RooStats::SPlot("sData","",data,&fitModel,RooArgList(Nsig,Nbkg));

  RooDataSet sigData("sigData","",&data,*data.get(),0,"Nsig_sw");
  RooDataSet bkgData("bkgData","",&data,*data.get(),0,"Nbkg_sw");
  */
  ws->import(data);
  //ws->import(sigData);
  //ws->import(bkgData);
  ws->import(fitModel);
  ws->import(*res);
  ws->import(*fmgg);
  ws->import(*fMR);

  return ws;
}

TString makeDoubleExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_dexp_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *alpha2 = new RooRealVar(tag+"_dexp_alpha2","#alpha_{2}",-0.1,-10,-0.0001);
  RooRealVar *f      = new RooRealVar(tag+"_dexp_f","f",0.3,0.0001,0.9999);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_dexp_Nbkg","N_{bkg}",10,1,1E9);

  RooExponential *exp1 = new RooExponential(tag+"_dexp_exp1","",mgg,*alpha1);
  RooExponential *exp2 = new RooExponential(tag+"_dexp_exp2","",mgg,*alpha2);

  RooAddPdf *add       = new RooAddPdf(tag+"_dexp","",*exp1,*exp2,*f);

  w.import(*(new RooExtendPdf(tag+"_dexp_ext","",*add,*Nbkg)));;

  return "dexp";
}

TString makeSingleExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_sexp_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_sexp_Nbkg","N_{bkg}",10,1,1E9);

  RooExponential *exp1 = new RooExponential(tag+"_sexp","",mgg,*alpha1);

  w.import(*(new RooExtendPdf(tag+"_sexp_ext","",*exp1,*Nbkg)));;

  return "sexp";
}

TString makeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_mexp_alpha1","#alpha_{1}",-1,-1000,-0.0000001);
  RooRealVar *m1 = new RooRealVar(tag+"_mexp_m1","m_{1}",1.,0.,10.);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_mexp_Nbkg","N_{bkg}",10,1,1E9);

  RooGenericPdf *mexp1 = new RooGenericPdf(tag+"_mexp","","exp(@0,*@1^@2)",RooArgList(*alpha1,mgg,*m1));

  w.import(*(new RooExtendPdf(tag+"_mexp_ext","",*mexp1,*Nbkg)));;

  return "mexp";
}

TString makeSinglePow(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_spow_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_spow_Nbkg","N_{bkg}",10,1,1E9);

  RooGenericPdf *pow1 = new RooGenericPdf(tag+"_spow","","@0^@1",RooArgList(mgg,*alpha1));

  w.import(*(new RooExtendPdf(tag+"_spow_ext","",*pow1,*Nbkg)));;

  return "spow";
}

TString makeDoublePow(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_dpow_alpha1","#alpha_{1}",1,0,10);
  RooRealVar *alpha2 = new RooRealVar(tag+"_dpow_alpha2","#alpha_{2}",2,0,10);
  RooRealVar *f      = new RooRealVar(tag+"_dpow_f","f",0.3,0.0001,0.9999);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_dpow_Nbkg","N_{bkg}",10,1,1E9);

  RooGenericPdf *pow1 = new RooGenericPdf(tag+"_dpow_pow1","","@0^@1",RooArgList(mgg,*alpha1));
  RooGenericPdf *pow2 = new RooGenericPdf(tag+"_dpow_pow2","","@0^@1",RooArgList(mgg,*alpha2));

  RooAddPdf *add       = new RooAddPdf(tag+"_dpow","",*pow1,*pow2,*f);

  w.import(*(new RooExtendPdf(tag+"_dpow_ext","",*add,*Nbkg)));

  return "dpow";
}




void scaleWithError(TH1D& hist,float scale) {
  for(int i=0;i<hist.GetNbinsX()+2;i++) {
    float v = hist.GetBinContent(i);

    hist.SetBinContent(i,v*scale);
    hist.SetBinError(i,sqrt(v)*scale);
  }
}


TCanvas * makeSigModelFit(TFile* gg,TFile* vbf, TFile* wz, TFile* tt, TString sel) {
  TH1D h_gg("h_gg","",40,110,170);
  TH1D h_vbf("h_vbf","",40,110,170);
  TH1D h_wz("h_wz","",40,110,170);
  TH1D h_tt("h_tt","",40,110,170);

  ((TTree*)gg->Get("SusyHggTree"))->Project("h_gg","mgg","("+sel+")");
  ((TTree*)vbf->Get("SusyHggTree"))->Project("h_vbf","mgg","("+sel+")");
  ((TTree*)wz->Get("SusyHggTree"))->Project("h_wz","mgg","("+sel+")");
  ((TTree*)tt->Get("SusyHggTree"))->Project("h_tt","mgg","("+sel+")");

  scaleWithError(h_gg,19.27*2.28E-3*19780/96290.);
  scaleWithError(h_vbf,1.578*2.28E-3*19780/99885.);
  scaleWithError(h_wz,0.7046*2.28E-3*19780/100320.);
  scaleWithError(h_tt,0.1293*2.28E-3*19780/93312.);

  TH1D* total = (TH1D*)h_gg.Clone("h_total");
  total->Add(&h_vbf);
  total->Add(&h_wz);
  total->Add(&h_tt);

  TF1* fit = new TF1("fit","gaus",120,130);
  fit->SetParameter(0,total->Integral());
  fit->SetParameter(1,125);
  fit->SetParameter(2,2);
  
  fit->SetParLimits(0,0.001,1E6);
  fit->SetParLimits(1,120,130);
  fit->SetParLimits(2,0.2,10);

  total->Fit(fit);
  TCanvas *c = new TCanvas("c","");

  total->SetXTitle("m_{#gamma#gamma} [GeV]");
  total->SetYTitle("Events / (1.5 GeV)");

  total->Draw("PE0");
  fit->Draw("SAME");

  return c;
  
}
