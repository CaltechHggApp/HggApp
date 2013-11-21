/*
  Methods to run the SUSY Hgg fits

*/

#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooCategory.h"

#include "TMath.h"
#include "TH2F.h"
#include "TString.h"

#include <vector>
#include <map>
#include "assert.h"

#define NUM_CPU 1

void RunHggFits(RooWorkspace * const ws,const char* tag);

TH2F* makeHistogram(const RooWorkspace* const ws, const char* dataName, const char*tag, const char* name);

void getScale(RooWorkspace* const ws,const char* tag);

std::vector<TString> defineCategories(RooWorkspace* const ws, const std::vector<TString>& mc_names);

TH2F* getSignalRegionHistogram(const TH2F& dataHist,const char* name);

void getMcNames(const RooWorkspace* const ws, std::vector<TString>& mc_names);

//entry function for all fitting
void RunFits(RooWorkspace * const ws,TFile *outputFile,const std::vector<TString>& mc_names) {
  std::vector<TString> cats = defineCategories(ws,mc_names);

  std::map<TString,float> xsec {
    {"gg_H_125",19.27},{"wz_H_125",0.7046},{"vbf_H_125",1.578},{"sms",0.3}
  };

  for(auto& cat : cats) {
    RunHggFits(ws,cat.Data());  
    TH2F* side = makeHistogram(ws,"data",cat.Data(),"Sideband");
    TH2F* signal = makeHistogram(ws,"data",cat.Data(),"Signal");
    side->SetName("data_"+cat+"_sideband");
    signal->SetName("data_"+cat+"_signal");

    getScale(ws,cat.Data());
    RooRealVar *scale = ws->var(Form("bkg_scale_factor_%s",cat.Data()));
    

    std::cout << "Scale: " << std::endl;
    scale->Print();

    TH2F* side_plus = (TH2F*)side->Clone("data_"+cat+"_sideband_up1sigma");
    TH2F* side_minus = (TH2F*)side->Clone("data_"+cat+"_sideband_down1sigma");
    
    side->Scale(scale->getVal());
    side_plus->Scale(scale->getVal()+scale->getError());
    side_minus->Scale(scale->getVal()-scale->getError());

    //deal with the signal MC
    TH2F* sm_tot = 0;
    for(auto& mc : mc_names) {
      TH2F* mc_hist = makeHistogram(ws,mc.Data(),cat.Data(),"All");
      mc_hist->SetName(mc+"_"+cat+"_signal");

      RooRealVar *n = ws->var("N_"+mc);

      assert(n->getVal()!=0);

      if(mc.First("sms")==-1) mc_hist->Scale(xsec[mc] * 2.28E-3 * 19780/n->getVal());
      else mc_hist->Scale(0.3 * 2.28E-3 * 19780/n->getVal());

      TH2F* signalRegion_sig = getSignalRegionHistogram(*mc_hist,mc+"_"+cat+"_signal_SigRegions");

      outputFile->cd();
      mc_hist->Write();
      signalRegion_sig->Write();
      if(mc.First("sms") != -1) continue;
      
      if(sm_tot==0) {
	sm_tot = (TH2F*)mc_hist->Clone("SMTot_"+cat+"_signal");		      
      }else{
	sm_tot->Add(mc_hist);
      }
    }
    TH2F* bkgSub = (TH2F*)signal->Clone("data_"+cat+"_signal_sidebandSub");
    TH2F* bkgSub1up = (TH2F*)signal->Clone("data_"+cat+"_signal_sidebandSub_up1sigma");
    TH2F* bkgSub1down = (TH2F*)signal->Clone("data_"+cat+"_signal_sidebandSub_down1sigma");
    RooRealVar *expected_sideband_integral = ws->var("sideband_fit_int_"+cat);

    std::cout << expected_sideband_integral->getVal()/side->Integral()*scale->getVal() << std::endl;
    std::cout << bkgSub->Integral() << std::endl;
    std::cout << expected_sideband_integral->getVal() << std::endl;
    std::cout << side->Integral() << std::endl;

    bkgSub->Add(side,-1 * expected_sideband_integral->getVal()/side->Integral()*scale->getVal());
    bkgSub1up->Add(side_plus,-1 * expected_sideband_integral->getVal()/side->Integral()*scale->getVal());
    bkgSub1down->Add(side_minus,-1 * expected_sideband_integral->getVal()/side->Integral()*scale->getVal());

    TH2F* sideRegion_data = getSignalRegionHistogram(*side,"data_"+cat+"_sideband_SigRegions");
    TH2F* signalRegion_data = getSignalRegionHistogram(*signal,"data_"+cat+"_signal_SigRegions");

    TH2F* signalRegion_data_sub = getSignalRegionHistogram(*bkgSub,"data_"+cat+"_signal_sidebandSub_SigRegions");
    TH2F* signalRegion_data_sub_1up = getSignalRegionHistogram(*bkgSub1up,"data_"+cat+"_signal_sidebandSub_up1sigma_SigRegions");
    TH2F* signalRegion_data_sub_1down = getSignalRegionHistogram(*bkgSub1down,"data_"+cat+"_signal_sidebandSub_down1sigma_SigRegions");

    TH2F* signalRegion_sig = getSignalRegionHistogram(*sm_tot,"SMTot_"+cat+"_signal_SigRegions");

    
    outputFile->cd();
    side->Write();
    side_plus->Write();
    side_minus->Write();
    signal->Write();
    sm_tot->Write();
    bkgSub->Write();
    bkgSub1up->Write();
    bkgSub1down->Write();

    signalRegion_data->Write();
    sideRegion_data->Write();

    signalRegion_data_sub->Write();
    signalRegion_data_sub_1up->Write();
    signalRegion_data_sub_1down->Write();
    signalRegion_sig->Write();

    //->Write();
    
  }
}



//run mgg fits using a double exponential
void RunHggFits(RooWorkspace * const ws,const char* tag) {
  RooRealVar *mgg = ws->var("mgg");

  //mgg->setRange(110,150);

  RooRealVar a1(Form("%s_alpha1",tag),"",-0.1,-1.,0.);
  RooRealVar a2(Form("%s_alpha2",tag),"",-0.1,-1.,0.);
  RooRealVar f (Form("%s_f",tag),"",0.1,0.,1.);

  RooExponential exp1(Form("%s_exp1",tag),"",*mgg,a1);
  RooExponential exp2(Form("%s_exp2",tag),"",*mgg,a2);
  RooAddPdf      dexp(Form("%s_dexp",tag),"",RooArgList(exp1,exp2),f);

  RooRealVar Nbkg(Form("%s_Nbkg",tag),"",1000,0,1e9);

  RooExtendPdf  bkgModel(Form("%s_bkgModel",tag),"",dexp,Nbkg);

  RooDataSet *ds = (RooDataSet*)ws->data(Form("data_%s",tag));
  bkgModel.fitTo(*ds,RooFit::Save(kFALSE),RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU));
  RooFitResult *res =   bkgModel.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU));
  res->SetName( Form("%s_fitResult",tag) );

  ws->import(bkgModel);
  ws->import(*res);
}

void getDataInRange(RooWorkspace * const ws,const char* tag,TString selection) {
  RooDataSet *data = (RooDataSet*)ws->data("data");

  data = (RooDataSet*)data->reduce(selection.Data());

  ws->import(*data,RooFit::Rename(tag));
}

TH2F* makeHistogram(const RooWorkspace* const ws, const char* dataName, const char* tag, const char* name) {
  RooRealVar* MR = ws->var("MR");
  RooRealVar* R  = ws->var("Rsq");

  MR->setRange(0,2000.);
  MR->setBins(50);

  R->setRange(0,1.0);
  R->setBins(25);

  TString fullname = Form("%s_%s_%s",dataName,tag,name);
  
  TString regionString = "(mgg > 123 && mgg < 127)";

  if( strcmp(name,"Sideband")==0) {
    regionString = "(mgg > 117 && mgg <121) || (mgg > 129 && mgg < 133)";
  }

  if(strncmp(name,"All",3)==0) regionString="1";

  return ((RooDataSet*)ws->data(Form("%s_%s",dataName,tag))->reduce(regionString))->createHistogram(*MR,*R,"",fullname);
}

void getScale(RooWorkspace* const ws,const char* tag) {

  RooRealVar *a1 = ws->var( Form("%s_alpha1",tag) );
  RooRealVar *a2 = ws->var( Form("%s_alpha2",tag) );
  RooRealVar *f = ws->var( Form("%s_f",tag) );
  RooRealVar *Nbkg = ws->var( Form("%s_Nbkg",tag) );

  

  RooRealVar *mgg = ws->var("mgg");

  assert(a1!=0 && a2!=0 && f!=0);

  RooRealVar min("min","",100);
  RooRealVar max("max","",180);

  RooRealVar min_side_low("min_side_low","",117);
  RooRealVar max_side_low("max_side_low","",121);
  RooRealVar min_side_high("min_side_high","",129);
  RooRealVar max_side_high("max_side_high","",133);
  RooRealVar min_signal("min_signal","",123);
  RooRealVar max_signal("max_signal","",127);
    
  RooFormulaVar integral_total("integral_total","","(@2/@0*(exp(@0*@4)-exp(@0*@3))+(1-@2)/@1*(exp(@1*@4)-exp(@1*@3)))",RooArgList(*a1,*a2,*f,min,max));
  RooFormulaVar integral_side_low("integral_side_low","","@5*(@2/@0*(exp(@0*@4)-exp(@0*@3))+(1-@2)/@1*(exp(@1*@4)-exp(@1*@3)))/@6",RooArgList(*a1,*a2,*f,min_side_low,max_side_low,*Nbkg,integral_total));
  RooFormulaVar integral_side_high("integral_side_high","","@5*(@2/@0*(exp(@0*@4)-exp(@0*@3))+(1-@2)/@1*(exp(@1*@4)-exp(@1*@3)))/@6",RooArgList(*a1,*a2,*f,min_side_high,max_side_high,*Nbkg,integral_total));
  RooFormulaVar integral_signal("integral_signal","","@5*(@2/@0*(exp(@0*@4)-exp(@0*@3))+(1-@2)/@1*(exp(@1*@4)-exp(@1*@3)))/@6",RooArgList(*a1,*a2,*f,min_signal,max_signal,*Nbkg,integral_total));

  RooFitResult *fitres = (RooFitResult*)ws->obj( Form("%s_fitResult",tag) );

  RooFormulaVar ratio("ratio","","@0/(@1+@2)",RooArgList(integral_signal,integral_side_low,integral_side_high));

  std::cout << Nbkg->GetName() << std::endl;
  std::cout << ratio.getVal() << std::endl << ratio.getPropagatedError(*fitres) << std::endl;

  RooRealVar bkg_scale(Form("bkg_scale_factor_%s",tag),"",ratio.getVal());
  bkg_scale.setError( ratio.getPropagatedError(*fitres) ); //deal with the correlations correctly

  RooRealVar sig_int(Form("signal_fit_int_%s",tag),"",integral_signal.getVal());
  sig_int.setError(integral_signal.getPropagatedError(*fitres));

  RooRealVar side_low_int(Form("sideband_low_fit_int_%s",tag),"",integral_side_low.getVal());
  side_low_int.setError(integral_side_low.getPropagatedError(*fitres));
  
  RooRealVar side_high_int(Form("sideband_high_fit_int_%s",tag),"",integral_side_high.getVal());
  side_high_int.setError(integral_side_high.getPropagatedError(*fitres));

  RooFormulaVar integral_side("integral_side","@0+@1",RooArgList(side_high_int,side_low_int));

  RooRealVar side_int(Form("sideband_fit_int_%s",tag),"",integral_side.getVal());
  side_int.setError(integral_side.getPropagatedError(*fitres));

  ws->import(bkg_scale);
  ws->import(ratio);

  ws->import(sig_int);
  ws->import(side_low_int);
  ws->import(side_high_int);
  ws->import(side_int);
}

void defineCategories(RooWorkspace* const ws, const char* dataName);

std::vector<TString> defineCategories(RooWorkspace* const ws,const std::vector<TString>& mc_names) {
  std::vector<TString> cats;
  cats.push_back("SingleMu");
  cats.push_back("SingleEle");
  cats.push_back("Hadronic_EBEB");
  //cats.push_back("Hadronic_NEBEB");

  defineCategories(ws,"data");
  for(auto& cat: mc_names) {
    defineCategories(ws,cat);
  }
  return cats;
}

void defineCategories(RooWorkspace* const ws, const char* dataName) {
  {

    RooDataSet *data = (RooDataSet*)ws->data(dataName);
    data = (RooDataSet*)data->reduce("mu1_pt > 20.");
    ws->import(*data,RooFit::Rename(Form("%s_SingleMu",dataName)) );
  }
  
  {
    RooDataSet *data = (RooDataSet*)ws->data(dataName);
    data = (RooDataSet*)data->reduce("mu1_pt <= 20. && ele1_pt > 20.");
    ws->import(*data,RooFit::Rename(Form("%s_SingleEle",dataName)) );
  }

  {
    RooDataSet *data = (RooDataSet*)ws->data(dataName);
    data = (RooDataSet*)data->reduce("mu1_pt <= 20. && ele1_pt <=20.");
    ws->import(*data,RooFit::Rename(Form("%s_Hadronic_EBEB",dataName)) );
  }
  /*  
  {
    RooDataSet *data = (RooDataSet*)ws->data(dataName);
    data = (RooDataSet*)data->reduce("mu1_pt <= 20. && ele1_pt <=20. && abs(pho1_eta) < 1.48 && abs(pho2_eta) < 1.48");
    ws->import(*data,RooFit::Rename(Form("%s_Hadronic_EBEB",dataName)) );
  }


  {
    RooDataSet *data = (RooDataSet*)ws->data(dataName);
    data = (RooDataSet*)data->reduce("mu1_pt <= 20. && ele1_pt <=20. && (abs(pho1_eta) > 1.48 || abs(pho2_eta) > 1.48)");
    ws->import(*data,RooFit::Rename(Form("%s_Hadronic_NEBEB",dataName)) );
  }
  */
  
      
}



TH2F* getSignalRegionHistogram(const TH2F& dataHist,const char* name) {

  //define the signal regions

  const int nXbins =5;
  double xBins[nXbins] = {0,200,400,1000,2000};
  const int nYbins =5;
  double yBins[nYbins] = {0,0.1,0.2,0.5,1.0};

  TH2F *signalHist = new TH2F(name,"",nXbins-1,xBins,nYbins-1,yBins);

  for(int iX=1;iX<dataHist.GetNbinsX();iX++) {
    for(int iY=1;iY<dataHist.GetNbinsY();iY++) {
      signalHist->Fill( dataHist.GetXaxis()->GetBinCenter(iX), dataHist.GetYaxis()->GetBinCenter(iY), dataHist.GetBinContent(iX,iY) );
    }
  }
  return signalHist;
}


void getMcNames(const RooWorkspace* const ws, std::vector<TString>& mc_names) {
  RooCategory *v = (RooCategory*)ws->obj("mc_names");
  int i=0;
  while(v->setIndex(i++)==0) {
    mc_names.push_back(v->getLabel());
  }
}
