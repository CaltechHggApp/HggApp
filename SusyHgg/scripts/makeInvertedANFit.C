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
#include "RooBernstein.h"
#include "RooMinuit.h"

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
#include "TMath.h"

#include <vector>
#include "assert.h"

#define DEBUG 1

TString makeDoubleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeSingleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeTripleExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeSinglePow(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makeDoublePow(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makePoly2(TString tag, RooRealVar& mgg,RooWorkspace& w);
TString makePoly3(TString tag, RooRealVar& mgg,RooWorkspace& w);

const int nTags=8;
TString alltags[nTags] = {"dexp","sexp","texp","mexp","spow","dpow","pol2","pol3"};

RooWorkspace* makeInvertedANFit(TTree* tree, float forceSigma=-1, bool constrainMu=false, float forceMu=-1) {
  RooWorkspace *ws = new RooWorkspace("ws","");

  std::vector< TString (*)(TString, RooRealVar&, RooWorkspace&) > bkgPdfList;
  bkgPdfList.push_back(makeDoubleExp);
#if DEBUG==0
  bkgPdfList.push_back(makeSingleExp);
  bkgPdfList.push_back(makeTripleExp);
  bkgPdfList.push_back(makeModExp);
  bkgPdfList.push_back(makeSinglePow);
  bkgPdfList.push_back(makeDoublePow);
  bkgPdfList.push_back(makePoly2);
#endif
  bkgPdfList.push_back(makePoly3);



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


  std::vector<TString> tags;
  //fit many different background models
  for(auto func = bkgPdfList.begin(); func != bkgPdfList.end(); func++) {
    TString tag = (*func)("bo",mgg,*ws);
    tags.push_back(tag);
    ws->pdf("bo_"+tag+"_model")->fitTo(data,RooFit::Strategy(0),RooFit::Extended(kTRUE));;
    RooFitResult* bres = ws->pdf("bo_"+tag+"_model")->fitTo(data,RooFit::Strategy(2),RooFit::Save(kTRUE),RooFit::Extended(kTRUE));
    bres->SetName("bo_"+tag+"_fitres");
    ws->import(*bres);
  }

  //fit all the S+B models
  for(auto func = bkgPdfList.begin(); func != bkgPdfList.end(); func++) {
    TString tag = (*func)("b",mgg,*ws);


    RooRealVar sigma("s_"+tag+"_sigma","",5,0,100);
    if(forceSigma!=-1) {
      sigma.setVal(forceSigma);
      sigma.setConstant(true);
    }
    RooRealVar mu("s_"+tag+"_mu","",126,120,130);
    if(forceMu!=-1) {
      mu.setVal(forceMu);
      mu.setConstant(true);
    }
    RooGaussian sig("s_"+tag+"_model","",mgg,mu,sigma);
    RooRealVar Nsig("sb_"+tag+"_Ns","",5,0,100);
    RooRealVar Nbkg("sb_"+tag+"_Nb","",100,0,100000);
    

    RooRealVar HiggsMass("HiggsMass","",125.7);
    RooRealVar HiggsMassError("HiggsMassError","",0.4);
    RooGaussian HiggsMassConstraint("HiggsMassConstraint","",mu,HiggsMass,HiggsMassError);


    RooAddPdf fitModel("sb_"+tag+"_model","",RooArgList( *ws->pdf("b_"+tag), sig ),RooArgList(Nbkg,Nsig));

    RooFitResult* sbres;
    if(constrainMu) {
      fitModel.fitTo(data,RooFit::Strategy(0),RooFit::Extended(kTRUE),RooFit::ExternalConstraints(RooArgSet(HiggsMassConstraint)));
      sbres = fitModel.fitTo(data,RooFit::Strategy(2),RooFit::Save(kTRUE),RooFit::Extended(kTRUE),RooFit::ExternalConstraints(RooArgSet(HiggsMassConstraint)));
    } else {
      fitModel.fitTo(data,RooFit::Strategy(0),RooFit::Extended(kTRUE));
      sbres = fitModel.fitTo(data,RooFit::Strategy(2),RooFit::Save(kTRUE),RooFit::Extended(kTRUE));
    }
    sbres->SetName("sb_"+tag+"_fitres");
    ws->import(*sbres);
    ws->import(fitModel);

    RooPlot *fmgg = mgg.frame();
    data.plotOn(fmgg);
    fitModel.plotOn(fmgg);
    ws->pdf("bo_"+tag+"_model")->plotOn(fmgg,RooFit::LineColor(kRed),RooFit::Range("Full"),RooFit::NormRange("Full"));
    fmgg->SetName(tag+"_frame");
    ws->import(*fmgg);
    delete fmgg;

    /*
    RooPlot *fNs = Nsig.frame(0,25);
    fNs->SetName(tag+"_Nsig_pll");
    RooAbsReal *nll = fitModel.createNLL(data,RooFit::NumCPU(4));
    RooMinuit(*nll).migrad();
    RooAbsReal *pll = nll->createProfile(Nsig);
    //nll->plotOn(fNs,RooFit::ShiftToZero(),RooFit::LineColor(kRed));
    pll->plotOn(fNs);
    ws->import(*fNs);
    delete fNs;
    */
  }

  RooRealVar minAIC("minAIC","",1E10);
  float aicExpSum=0;
  //compute AIC stuff
  const int nModels=2;
  TString models[nModels] = {"bo","sb"};
  for(auto t = tags.begin(); t!=tags.end(); t++) {
    for(int iModel=0; iModel<nModels; iModel++) {
      RooAbsPdf *p = ws->pdf(models[iModel]+"_"+*t+"_model");
      RooFitResult *res = (RooFitResult*)ws->obj(models[iModel]+"_"+*t+"_fitres");
      assert(p!=0);
      assert(res!=0);
      RooRealVar k(models[iModel]+"_"+*t+"_k","",p->getParameters(RooArgSet(mgg))->getSize()-(iModel==1?2:0));
      RooRealVar nll(models[iModel]+"_"+*t+"_minNll","",res->minNll());
      RooFormulaVar AIC(models[iModel]+"_"+*t+"_AIC","2*@0+2*@1",RooArgSet(nll,k));
      ws->import(AIC);
      if(AIC.getVal() < minAIC.getVal()) {
	minAIC.setVal(AIC.getVal());
      }
      aicExpSum+=TMath::Exp(-0.5*AIC.getVal()); //we will need this precomputed  for the next step
    }
  }
  ws->import(minAIC);
  //compute the AIC weight
  for(auto t = tags.begin(); t!=tags.end(); t++) {
    cout << *t << std::endl;
    for(int iModel=0; iModel<nModels; iModel++) {
      RooFormulaVar *AIC = (RooFormulaVar*)ws->obj(models[iModel]+"_"+*t+"_AIC");
      RooRealVar AICw(models[iModel]+"_"+*t+"_AICWeight","",TMath::Exp(-0.5*(AIC->getVal()-minAIC.getVal()))/aicExpSum/TMath::Exp(0.5*minAIC.getVal()));
      RooRealVar AICp(models[iModel]+"_"+*t+"_AICProb","",TMath::Exp(-0.5*(AIC->getVal()-minAIC.getVal())));
      ws->import(AICw);
      std::cout << "\t" << models[iModel] << " &  " << AIC->getVal()-minAIC.getVal() << "  &  " << AICw.getVal()  << "\\\\" << std::endl;
      
    }
  }
  

  return ws;
}

float getDLL(RooWorkspace* w, TString tag) {
  //RooFitResult::sexp_b_fitres

  RooFitResult *b = (RooFitResult*)w->obj(tag+"_b_fitres");
  RooFitResult *sb = (RooFitResult*)w->obj(tag+"_sb_fitres");

  return -2*(sb->minNll()-b->minNll());
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

  w.import(*(new RooExtendPdf(tag+"_dexp_model","",*add,*Nbkg)));;

  return "dexp";
}

TString makeSingleExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_sexp_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_sexp_Nbkg","N_{bkg}",10,1,1E9);

  RooExponential *exp1 = new RooExponential(tag+"_sexp_model","",mgg,*alpha1);

  w.import(*(new RooExtendPdf(tag+"_sexp_model","",*exp1,*Nbkg)));;

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

  w.import(*(new RooExtendPdf(tag+"_texp_model","",*add,*Nbkg)));;

  return "texp";
}

TString makeModExp(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_mexp_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *m1 = new RooRealVar(tag+"_mexp_m1","m_{1}",1.,0.,10.);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_mexp_Nbkg","N_{bkg}",10,1,1E9);

  RooGenericPdf *mexp1 = new RooGenericPdf(tag+"_mexp","","exp(@0,*@1^@2)",RooArgList(*alpha1,mgg,*m1));

  w.import(*(new RooExtendPdf(tag+"_mexp_model","",*mexp1,*Nbkg)));;

  return "mexp";
}

TString makeSinglePow(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *alpha1 = new RooRealVar(tag+"_spow_alpha1","#alpha_{1}",-1,-10,-0.0001);
  RooRealVar *Nbkg   = new RooRealVar(tag+"_spow_Nbkg","N_{bkg}",10,1,1E9);

  RooGenericPdf *pow1 = new RooGenericPdf(tag+"_spow","","@0^@1",RooArgList(mgg,*alpha1));

  w.import(*(new RooExtendPdf(tag+"_spow_model","",*pow1,*Nbkg)));;

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

  w.import(*(new RooExtendPdf(tag+"_dpow_model","",*add,*Nbkg)));

  return "dpow";
}

TString makePoly2(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *pC = new RooRealVar(tag+"_pol2_pC","C",1,-10,10);
  RooRealVar *p0 = new RooRealVar(tag+"_pol2_p0","p_0",0,-10,10);
  RooRealVar *p1 = new RooRealVar(tag+"_pol2_p1","p_1",0,-10,10);

  RooRealVar *Nbkg   = new RooRealVar(tag+"_pol2_Nbkg","N_{bkg}",10,1,1E9);

  RooFormulaVar *pCmod = new RooFormulaVar(tag+"_pol2_pCmod","@0*@0",*pC);
  RooFormulaVar *p0mod = new RooFormulaVar(tag+"_pol2_p0mod","@0*@0",*p0);
  RooFormulaVar *p1mod = new RooFormulaVar(tag+"_pol2_p1mod","@0*@0",*p1);

  RooBernstein* bern = new RooBernstein(tag+"_pol2","",mgg,RooArgList(*pCmod,*p0mod,*p1mod));

  w.import(*(new RooExtendPdf(tag+"_pol2_model","",*bern,*Nbkg)));

  return "pol2";
}

TString makePoly3(TString tag, RooRealVar& mgg,RooWorkspace& w) {
  RooRealVar *pC = new RooRealVar(tag+"_pol3_pC","C",1,-10,10);
  RooRealVar *p0 = new RooRealVar(tag+"_pol3_p0","p_0",0,-10,10);
  RooRealVar *p1 = new RooRealVar(tag+"_pol3_p1","p_1",0,-10,10);
  RooRealVar *p2 = new RooRealVar(tag+"_pol3_p2","p_2",0,-10,10);

  RooRealVar *Nbkg   = new RooRealVar(tag+"_pol3_Nbkg","N_{bkg}",10,1,1E9);

  RooFormulaVar *pCmod = new RooFormulaVar(tag+"_pol3_pCmod","@0*@0",*pC);
  RooFormulaVar *p0mod = new RooFormulaVar(tag+"_pol3_p0mod","@0*@0",*p0);
  RooFormulaVar *p1mod = new RooFormulaVar(tag+"_pol3_p1mod","@0*@0",*p1);
  RooFormulaVar *p2mod = new RooFormulaVar(tag+"_pol3_p2mod","@0*@0",*p2);

  RooBernstein* bern = new RooBernstein(tag+"_pol3","",mgg,RooArgList(*pCmod,*p0mod,*p1mod,*p2mod));

  w.import(*(new RooExtendPdf(tag+"_pol3_model","",*bern,*Nbkg)));

  return "pol3";
}
