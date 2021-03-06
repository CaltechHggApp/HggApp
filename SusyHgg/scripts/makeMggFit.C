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
#include "RooFormulaVar.h"

#include "TLatex.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLegend.h"

TString makeDoubleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeSingleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeTripleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeSinglePow(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeDoublePow(TString tag, RooRealVar& mgg,RooWorkspace& w);


RooWorkspace* makeMggFit(TTree* tree,bool doBkg=true,TString plotBkg="", float forceSigma=-1, float forceNsig=-1, float forceMu=-1,float forceA1=1,float forceA2=1, float forceF=-1) {
  RooWorkspace *ws = new RooWorkspace("ws","");

  std::vector< TString (*)(TString, RooRealVar&, RooWorkspace&) > bkgPdfList;
  bkgPdfList.push_back(makeDoubleExp);
  bkgPdfList.push_back(makeSingleExp);
  bkgPdfList.push_back(makeTripleExp);
  bkgPdfList.push_back(makeModExp);
  bkgPdfList.push_back(makeSinglePow);
  bkgPdfList.push_back(makeDoublePow);



  RooRealVar mgg("mgg","m_{#gamma#gamma}",103,160,"GeV");
  mgg.setBins(38);

  mgg.setRange("sideband_low", 103,120);
  mgg.setRange("sideband_high",131,160);
  mgg.setRange("signal",120,131);

  RooRealVar MR("MR","",0,3000,"GeV");
  MR.setBins(60);
  
  RooRealVar Rsq("t1Rsq","",0,1,"GeV");
  Rsq.setBins(20);

  RooRealVar hem1_M("hem1_M","",-1,2000,"GeV");
  hem1_M.setBins(40);

  RooRealVar hem2_M("hem2_M","",-1,2000,"GeV");
  hem2_M.setBins(40);

  RooRealVar ptgg("ptgg","p_{T}^{#gamma#gamma}",0,500,"GeV");
  ptgg.setBins(50);

  RooDataSet data("data","",tree,RooArgSet(mgg,MR,Rsq,hem1_M,hem2_M,ptgg));

  RooRealVar alpha1("alpha1","",-1,-10,-0.001);
  if(forceA1!=1) {
    alpha1.setVal(forceA1);
    alpha1.setConstant(kTRUE);
  }
  RooRealVar alpha2("alpha2","",-1,-10,-0.001);
  if(forceA2!=1) {
    alpha2.setVal(forceA2);
    alpha2.setConstant(kTRUE);
  }
  RooRealVar f("f","",0.1,0.001,0.999);
  if(forceF!=-1) {
    f.setVal(forceA1);
    f.setConstant(kTRUE);
  }
  RooExponential exp1("exp1","",mgg,alpha1);
  RooExponential exp2("exp2","",mgg,alpha2);
  
  RooAddPdf bkgModel("bkgModel","",exp1,exp2,f);

  RooRealVar mu("mu","",125,120,130);
  if(forceMu>0) {
    mu.setRange(110,170);
    mu.setVal(forceMu);
    mu.setConstant(true);
  }
  RooRealVar sigma("sigma","",1,0.2,6);
  if(forceSigma>0) {
    sigma.setVal(forceSigma);
    sigma.setConstant(true);
  }

  RooRealVar Nbkg("Nbkg","",10,1,1E9);
  RooRealVar Nsig("Nsig","",1,0,1E4);  

  if(forceNsig!=-1) {
    Nsig.setVal(forceNsig);
    Nsig.setConstant(true);
  }

  RooGaussian sigModel("sigModel","",mgg,mu,sigma);

  RooAddPdf fitModel("fitModel","",RooArgList(bkgModel,sigModel),RooArgList(Nbkg,Nsig));

  fitModel.fitTo(data,RooFit::Strategy(0));
  RooFitResult *res = fitModel.fitTo(data,RooFit::Strategy(2),RooFit::Save(true),RooFit::Extended(kTRUE));
  res->SetName("fitres");

  //fit many different background models
  if(doBkg) {
    for(auto func = bkgPdfList.begin(); func != bkgPdfList.end(); func++) {
      TString tag = (*func)("bkgPdf",mgg,*ws);
      ws->pdf("bkgPdf_"+tag+"_ext")->fitTo(data,RooFit::Strategy(0),RooFit::Extended(kTRUE));
      RooFitResult* bres = ws->pdf("bkgPdf_"+tag+"_ext")->fitTo(data,RooFit::Strategy(2),RooFit::Save(kTRUE),RooFit::Extended(kTRUE));
      bres->SetName("fitres_bkgPdf_"+tag);
      ws->import(*bres);
    }
  }

  RooAbsReal *bkgInt = bkgModel.createIntegral(mgg,RooFit::NormSet(mgg),RooFit::Range("signal"));
  RooFormulaVar bkgTot("bkgTot","@0*@1",RooArgSet(*bkgInt,Nbkg));


  RooPlot *fmgg = mgg.frame();
  fmgg->SetName("mgg_frame");
  data.plotOn(fmgg);
  fitModel.plotOn(fmgg);

  TLatex muL(0.62,0.8,Form("#mu = %0.1f GeV",mu.getVal()));
  muL.SetTextSize(0.042);
  muL.SetNDC();

  TLatex sigL(0.62,0.72,Form("#sigma = %0.1f GeV",sigma.getVal()));
  sigL.SetTextSize(0.042);
  sigL.SetNDC();

  TLatex nL(0.62,0.64,Form("N_{sig} = %0.1f",Nsig.getVal()));
  nL.SetTextSize(0.042);
  nL.SetNDC();

  TLatex nB(0.62,0.56,Form("N_{bkg} = %0.2f",bkgTot.getVal()));
  nB.SetTextSize(0.042);
  nB.SetNDC();

  TLatex prelim(0.14,0.96,"CMS Preliminary");
  prelim.SetTextSize(0.042);
  prelim.SetNDC();
  
  TLatex dataL(0.6,0.96,"#sqrt{s} = 8 TeV  L = 19.8 fb^{-1}");
  dataL.SetTextSize(0.042);
  dataL.SetNDC();

  fmgg->addObject(&muL);
  fmgg->addObject(&sigL);
  fmgg->addObject(&nL);
  fmgg->addObject(&nB);

  fmgg->addObject(&prelim);
  fmgg->addObject(&dataL);


  if(plotBkg!="" && doBkg) {
    ws->pdf(Form("bkgPdf_%s_ext",plotBkg.Data()))->plotOn(fmgg,RooFit::LineColor(kRed));
    TH1D *bkg_dummy = new TH1D("bkg_dummy","",2,0,1);
    TH1D *sig_dummy = new TH1D("sig_dummy","",2,0,1);
    bkg_dummy->SetLineColor(kRed);
    sig_dummy->SetLineColor(kBlue);
    bkg_dummy->SetLineWidth(2);
    sig_dummy->SetLineWidth(2);
    TLegend *leg = new TLegend(0.62,0.85,0.85,0.95);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(sig_dummy,"S+B","l");
    leg->AddEntry(bkg_dummy,"B only","l");
    fmgg->addObject(leg);
    std::cout << "\n\n\nHypothesis Testing: S+B vs B only:" << std::endl;
    std::cout << "sqrt(2\Delta LL) =  " << sqrt( 2*( ((RooFitResult*)ws->obj(Form("fitres_bkgPdf_%s",plotBkg.Data())))->minNll() - res->minNll() ) ) << std::endl;
    std::cout << "\n\n\n" << std::endl;
  }
  
  alpha1.setConstant();
  alpha2.setConstant();
  sigma.setConstant();
  mu.setConstant();

  RooStats::SPlot* sdata = new RooStats::SPlot("sData","",data,&fitModel,RooArgList(Nsig,Nbkg));

  RooDataSet sigData("sigData","",&data,*data.get(),0,"Nsig_sw");
  RooDataSet bkgData("bkgData","",&data,*data.get(),0,"Nbkg_sw");
  
  std::cout << "Nbkg: " << Nbkg.getVal() << " \pm " << Nbkg.getError() << std::endl;
  std::cout << "Nbkg: " << bkgInt->getVal() << " \pm " << bkgInt->getPropagatedError(*res) << std::endl;

  std::cout << "BKG ERROR:  " << bkgTot.getPropagatedError(*res) << std::endl;

  ws->import(data);
  ws->import(sigData);
  ws->import(bkgData);
  ws->import(fitModel);
  ws->import(*res);
  ws->import(*fmgg);

  return ws;
}

TString makeDoubleExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_dexp_alpha1","#alpha_{1}",-0.03,-4.,4.);
  RooRealVar *alpha2 = new RooRealVar(tag+"_dexp_alpha2","#alpha_{2}",-0.01,-4.,4.);

  RooFormulaVar *asq1 = new RooFormulaVar(tag+"_dexp_asq1","","-1*@0^2",*alpha1);
  RooFormulaVar *asq2 = new RooFormulaVar(tag+"_dexp_asq2","","-1*@0^2",*alpha2);

  RooRealVar *f      = new RooRealVar(tag+"_dexp_f","f",0.3,0.001,0.999);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_dexp_Nbkg","N_{bkg}",10,0.001,1E9);

  //RooFormulaVar *fn = new RooFormulaVar(tag+"_dexp_fn","","0.5*(tanh(@0)+1)",*f);

  //RooExponential *exp1 = new RooExponential(tag+"_dexp_exp1","",mgg,*alpha1);
  //RooExponential *exp2 = new RooExponential(tag+"_dexp_exp2","",mgg,*alpha2);

  RooExponential *exp1 = new RooExponential(tag+"_dexp_exp1","",mgg,*asq1);
  RooExponential *exp2 = new RooExponential(tag+"_dexp_exp2","",mgg,*asq2);

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

TString makeTripleExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_texp_alpha1","#alpha_{1}",-0.03,-4.,4.);
  RooRealVar *alpha2 = new RooRealVar(tag+"_texp_alpha2","#alpha_{2}",-0.01,-4.,4.);
  RooRealVar *alpha3 = new RooRealVar(tag+"_texp_alpha3","#alpha_{3}",-0.0001,-4.,4.);

  RooFormulaVar *asq1 = new RooFormulaVar(tag+"_texp_asq1","","-1*@0^2",*alpha1);
  RooFormulaVar *asq2 = new RooFormulaVar(tag+"_texp_asq2","","-1*@0^2",*alpha2);
  RooFormulaVar *asq3 = new RooFormulaVar(tag+"_texp_asq2","","-1*@0^2",*alpha2);

  RooRealVar *f1      = new RooRealVar(tag+"_texp_f1","f1",0.3,0.001,0.999);
  RooRealVar *f2      = new RooRealVar(tag+"_texp_f2","f2",0.01,0.001,0.999);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_texp_Nbkg","N_{bkg}",10,0.001,1E9);

  RooExponential *exp1 = new RooExponential(tag+"_texp_exp1","",mgg,*asq1);
  RooExponential *exp2 = new RooExponential(tag+"_texp_exp2","",mgg,*asq2);
  RooExponential *exp3 = new RooExponential(tag+"_texp_exp3","",mgg,*asq3);

  RooAddPdf *add       = new RooAddPdf(tag+"_texp","",RooArgSet(*exp1,*exp2,*exp3),RooArgSet(*f1,*f2));

  w.import(*(new RooExtendPdf(tag+"_texp_ext","",*add,*Nbkg)));;

  return "texp";
}

TString makeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_mexp_alpha1","#alpha_{1}",-1,-10,-0.0001);
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


TCanvas * makeSigModelFit(TFile* gg,TFile* vbf, TFile* wz, TFile* tt, TString sel,TF1* fit) {
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

  fit = new TF1("smhfit","gaus",120,130);
  fit->SetParameter(0,total->Integral());
  fit->SetParameter(1,125);
  fit->SetParameter(2,2);
  
  fit->SetParLimits(0,0.00001,1E6);
  fit->SetParLimits(1,120,130);
  fit->SetParLimits(2,0.2,10);

  total->Fit(fit);
  TCanvas *c = new TCanvas("c","");

  total->SetXTitle("m_{#gamma#gamma} [GeV]");
  total->SetYTitle("Events / (1.5 GeV)");

  total->Draw("PE0");
  //fit->Draw("SAME");

  return c;
}

RooRealVar* MakeFit(RooDataSet* data,RooRealVar& mgg,float forceMu,float forceSigma) {

  RooRealVar alpha1("alpha1","",-1,-10,-0.001);
  RooRealVar alpha2("alpha2","",-1,-10,-0.001);
  RooRealVar f("f","",0.1,0.001,0.999);
  RooExponential exp1("exp1","",mgg,alpha1);
  RooExponential exp2("exp2","",mgg,alpha2);
  
  RooAddPdf bkgModel("bkgModel","",exp1,exp2,f);


  RooRealVar mu("mu","",125,120,130);
  if(forceMu>0) {
    mu.setRange(100,170);
    mu.setVal(forceMu);
    mu.setConstant(true);
  }
  RooRealVar sigma("sigma","",1,0.2,6);
  if(forceSigma>0) {
    sigma.setVal(forceSigma);
    sigma.setConstant(true);
  }

  RooRealVar Nbkg("Nbkg","",10,1,1E5);
  RooRealVar *Nsig = new RooRealVar("Nsig","",1,0,1E4);

  RooGaussian sigModel("sigModel","",mgg,mu,sigma);

  RooAddPdf fitModel("fitModel","",RooArgList(bkgModel,sigModel),RooArgList(Nbkg,*Nsig));

  fitModel.fitTo(*data,RooFit::Strategy(0),RooFit::PrintLevel(-1),RooFit::Warnings(-1),RooFit::Verbose(0),RooFit::PrintEvalErrors(-1));
  fitModel.fitTo(*data,RooFit::Strategy(2),RooFit::PrintLevel(-1),RooFit::Warnings(-1),RooFit::Verbose(0),RooFit::PrintEvalErrors(-1));
  

  return Nsig;
}


float getPV(TRandom3& rng, RooRealVar& Nbkg, RooRealVar& mgg, RooAbsPdf& bkgModel,float obs, float obsE,float forceMu,float forceSigma) {
  int N=10;
  int count=0;

  while(count < 10) {
    count =0;
    N*=10;
    for(int i=0;i<N;i++) {
      int Ntoy = rng.Gaus(Nbkg.getVal(),Nbkg.getError());
      RooDataSet* toyData = bkgModel.generate(mgg,Ntoy);
      
      RooRealVar* fitted = MakeFit(toyData,mgg,forceMu,forceSigma);

      float Nsig = rng.Gaus(obs,obsE);
      float Nfitsig = rng.Gaus(fitted->getVal(),fitted->getError());

      if(TMath::Abs(Nfitsig) >= TMath::Abs(Nsig)) count+=1;
      delete toyData;
      delete fitted;
    }
  }
  return float(count)/N;
}


TGraph* makePVScan(TTree* tree, float sigma, int mggMin,int mggMax,TGraphErrors* obsG,TGraphErrors* bkgG,
		   float forceA1=1,float forceA2=1, float forceF=-1) {
  TGraph* g = new TGraph(mggMax-mggMin+1);

  obsG = new TGraphErrors(mggMax-mggMin+1);
  bkgG = new TGraphErrors(mggMax-mggMin+1);

  TRandom3 rng(0);
  for(int mgg = mggMin; mgg <= mggMax; mgg++) {
    RooWorkspace* fit = makeMggFit(tree,sigma*mgg/125.,mgg,forceA1,forceA2,forceF);
    fit->Print();
    float obs = fit->var("Nsig")->getVal();
    float obsE = fit->var("Nsig")->getError();
    //std::cout << obs << " +- " << obsE << std::endl;

    RooAbsPdf* bkgpdf = fit->pdf("bkgModelOnly");
    //RooFitResult* fitr = (RooFitResult*)fit->obj("fitres");
    RooRealVar* rrvmgg = fit->var("mgg");
    //rrvmgg->setRange("sig",125-2*sigma,125+2*sigma);
    //RooAbsReal* bkgintegral = bkgpdf->createIntegral(*rrvmgg,RooFit::NormSet(*rrvmgg),RooFit::Range("sig"));
    RooRealVar* bkg = fit->var("NbkgOnly");
    float bkgE = 0; //integral->getPropagatedError(*fitr);

    //RooAbsReal* sigintegral = fitpdf->createIntegral(*rrvmgg,RooFit::NormSet(*rrvmgg),RooFit::Range("sig"));
    //float obs = sigintegral->getVal() * (fit->var("Nbkg")->getVal()+fit->var("Nsig")->getVal());
    //float obsE=0;

    //std::cout << bkg << " +- " << bkgE << std::endl;
    float pv = getPV(rng,*bkg,*rrvmgg,*bkgpdf,obs,obsE,mgg,sigma);

    g->SetPoint(mgg-mggMin,mgg,pv);
    
    obsG->SetPoint(mgg-mggMin,mgg,obs);
    obsG->SetPointError(mgg-mggMin,0,obsE);

    bkgG->SetPoint(mgg-mggMin,mgg,bkg->getVal());
    bkgG->SetPointError(mgg-mggMin,0,bkg->getError());

    delete fit;
  }
  return g;
}
