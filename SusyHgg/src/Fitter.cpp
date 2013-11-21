#include "Fitter.hpp"

#include "TObjArray.h"
#include "TMath.h"

#define NUM_CPU 1

Fitter::Fitter(TString inputFileName,TString outputFileName) {
  inputFile = new TFile(inputFileName);
  ws = (RooWorkspace*)inputFile->Get("susy_hgg_workspace");

  outputFile = new TFile(outputFileName,"RECREATE");

  TObjArray *mc_names_array = (TObjArray*)inputFile->Get("mc_name_array");
  for(int i=0;i<mc_names_array->GetEntries();i++) {
    this->addMcName( ((TObjString*)mc_names_array->At(i))->String() );
  }

  //
  defineCats();
}

Fitter::~Fitter() {
  inputFile->Close();
  outputFile->Close();
  delete ws;

}

void Fitter::Run() {
  for(int catIndex=0;catIndex < cats.size();catIndex++) {
    buildHistograms(catIndex);
    runFits(catIndex);
    computeScaleFactor(catIndex);
    buildSubtractedHistograms(catIndex);
  }
  saveAll();
}

void Fitter::buildHistograms(int catIndex) {
  std::cout << "buildHistograms " << catIndex << std::endl;

  save_histograms[Form("data_%s_Signal",cats.at(catIndex).Data())] = build_histogram("data",catIndex,kSignal);
  save_histograms[Form("data_%s_Sideband",cats.at(catIndex).Data())] = build_histogram("data",catIndex,kSideband);
  save_histograms[Form("data_%s_All",cats.at(catIndex).Data())] = build_histogram("data",catIndex,kAll);
  
  TH2F* SMTot=0;
  for(auto& mc_name : sm_higgs_names) {
    RooRealVar *n = ws->var("N_"+mc_name);
    save_histograms[Form("%s_%s_Signal",mc_name.Data(),cats.at(catIndex).Data())] = build_histogram(mc_name,catIndex,kSignal,xsecs[mc_name] * 2.28E-3 * 19780/n->getVal());
    if(SMTot==0) SMTot = (TH2F*) save_histograms[Form("%s_%s_Signal",mc_name.Data(),cats.at(catIndex).Data())]->Clone("SMTot_"+cats.at(catIndex)+"_Signal");
    else SMTot->Add(save_histograms[Form("%s_%s_Signal",mc_name.Data(),cats.at(catIndex).Data())]);
  }
  save_histograms[Form("SMTot_%s_Signal",cats.at(catIndex).Data())] = SMTot;
  
  for(auto& mc_name : sms_names) {
    RooRealVar *n = ws->var("N_"+mc_name);
    save_histograms[Form("%s_%s_Signal",mc_name.Data(),cats.at(catIndex).Data())] = build_histogram(mc_name,catIndex,kSignal,1.0 * 2.28E-3 * 19780/n->getVal());      
  }
  
 }


void Fitter::runFits(int catIndex) {
  std::cout << "runFits " << catIndex << std::endl;

  RooRealVar *mgg = ws->var("mgg");

  TString tag = cats.at(catIndex);

  //mgg->setRange(110,150);

  RooRealVar a1(Form("%s_alpha1",tag.Data()),"",-0.1,-1.,0.);
  RooRealVar a2(Form("%s_alpha2",tag.Data()),"",-0.1,-1.,0.);
  RooRealVar f (Form("%s_f",tag.Data()),"",0.1,0.,1.);

  RooExponential exp1(Form("%s_exp1",tag.Data()),"",*mgg,a1);
  RooExponential exp2(Form("%s_exp2",tag.Data()),"",*mgg,a2);
  RooAddPdf      dexp(Form("%s_dexp",tag.Data()),"",RooArgList(exp1,exp2),f);

  RooRealVar Nbkg(Form("%s_Nbkg",tag.Data()),"",1000,0,1e9);

  RooExtendPdf  bkgModel(Form("%s_bkgModel",tag.Data()),"",dexp,Nbkg);

  RooDataSet *ds = (RooDataSet*)ws->data(Form("data_%s",tag.Data()));
  bkgModel.fitTo(*ds,RooFit::Save(kFALSE),RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU));
  RooFitResult *res =   bkgModel.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU));
  res->SetName( Form("%s_fitResult",tag.Data()) );

  ws->import(bkgModel);
  ws->import(*res);
}

void Fitter::computeScaleFactor(int catIndex) {
  std::cout << "computeScaleFactor " << catIndex << std::endl;
  TString tag = cats.at(catIndex);
  RooRealVar *a1 = ws->var( Form("%s_alpha1",tag.Data()) );
  RooRealVar *a2 = ws->var( Form("%s_alpha2",tag.Data()) );
  RooRealVar *f = ws->var( Form("%s_f",tag.Data()) );
  RooRealVar *Nbkg = ws->var( Form("%s_Nbkg",tag.Data()) );

  

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

  RooFitResult *fitres = (RooFitResult*)ws->obj( Form("%s_fitResult",tag.Data()) );

  RooFormulaVar ratio("ratio","","@0/(@1+@2)",RooArgList(integral_signal,integral_side_low,integral_side_high));

  std::cout << Nbkg->GetName() << std::endl;
  std::cout << ratio.getVal() << std::endl << ratio.getPropagatedError(*fitres) << std::endl;

  RooRealVar bkg_scale(Form("bkg_scale_factor_%s",tag.Data()),"",ratio.getVal());
  bkg_scale.setError( ratio.getPropagatedError(*fitres) ); //deal with the correlations correctly

  RooRealVar sig_int(Form("signal_fit_int_%s",tag.Data()),"",integral_signal.getVal());
  sig_int.setError(integral_signal.getPropagatedError(*fitres));

  RooRealVar side_low_int(Form("sideband_low_fit_int_%s",tag.Data()),"",integral_side_low.getVal());
  side_low_int.setError(integral_side_low.getPropagatedError(*fitres));
  
  RooRealVar side_high_int(Form("sideband_high_fit_int_%s",tag.Data()),"",integral_side_high.getVal());
  side_high_int.setError(integral_side_high.getPropagatedError(*fitres));

  RooFormulaVar integral_side("integral_side","@0+@1",RooArgList(side_high_int,side_low_int));

  RooRealVar side_int(Form("sideband_fit_int_%s",tag.Data()),"",integral_side.getVal());
  side_int.setError(integral_side.getPropagatedError(*fitres));

  ws->import(bkg_scale);
  ws->import(ratio);

  ws->import(sig_int);
  ws->import(side_low_int);
  ws->import(side_high_int);
  ws->import(side_int);

  sideband_integrals.at(catIndex) = std::make_pair(side_int.getVal(),side_int.getError());
  signal_integrals.at(catIndex)   = std::make_pair(sig_int.getVal(), sig_int.getError());
  sideband_to_signal_scale_factors.at(catIndex) = std::make_pair(bkg_scale.getVal(),bkg_scale.getError());
}

void Fitter::buildSubtractedHistograms(int catIndex) {
  std::cout << "buildSubtractedHistogram " << catIndex << std::endl;
  TH2F* signal   = save_histograms[ "data_"+cats.at(catIndex)+"_Signal"];
  TH2F* sideband = save_histograms[ "data_"+cats.at(catIndex)+"_Sideband"];

  TH2F* bkg_sub = (TH2F*)signal->Clone("data_"+cats.at(catIndex)+"_Signal_sidebandSub");
  //bkg_sub->Add(sideband,-1*sideband_to_signal_scale_factors.at(catIndex).first);
  bkg_sub->Add(sideband,-1*signal_integrals.at(catIndex).first/sideband->Integral());

  for(int iX=1;iX < bkg_sub->GetNbinsX();iX++) {
    for(int iY=1;iY < bkg_sub->GetNbinsY();iY++) {
      bkg_sub->SetBinError(iX,iY,TMath::Sqrt(signal->GetBinContent(iX,iY)+sideband->GetBinContent(iX,iY)+1));
    }
  }

  save_histograms["data_"+cats.at(catIndex)+"_Signal_sidebandSub"] = bkg_sub;
}

TH2F* Fitter::getSignalRegionHistogram(const TH2F& hist,const char* name) {
  const int nXbins =5;
  double xBins[nXbins] = {0,200,400,1000,2000};
  const int nYbins =5;
  double yBins[nYbins] = {0,0.1,0.2,0.5,1.0};

  TH2F *signalHist = new TH2F(name,"",nXbins-1,xBins,nYbins-1,yBins);

  for(int iX=1;iX<hist.GetNbinsX();iX++) {
    for(int iY=1;iY<hist.GetNbinsY();iY++) {
      signalHist->Fill( hist.GetXaxis()->GetBinCenter(iX), hist.GetYaxis()->GetBinCenter(iY), hist.GetBinContent(iX,iY) );
    }
  }

  RooRealVar sigregion("signal_region","",-0.5,nXbins*nYbins-0.5);
  sigregion.setBins(nXbins*nYbins);

  RooDataHist rdh( Form("%s_hist",name),"",sigregion);

  for(int iX=0;iX<nXbins;iX++) {
    for(int iY=0;iY<nYbins;iY++) {
    
      sigregion.setVal( iX*4+iY );
      rdh.set(sigregion,signalHist->GetBinContent(iX+1,iY+1));
    }
  }

  RooHistPdf rhp( Form("%s_pdf",name), "", sigregion, rdh);


  ws->import(sigregion);
  ws->import(rhp);

  return signalHist;

}

void Fitter::addMcName(TString name) {
  if(name.Index("_H_")!=-1) sm_higgs_names.push_back(name);
  else sms_names.push_back(name);
}


TH2F* Fitter::build_histogram(TString name,int catIndex,region reg,float weight) {
  RooRealVar* MR = ws->var("MR");
  RooRealVar* R  = ws->var("Rsq");

  MR->setRange(0,2000.);
  MR->setBins(50);

  R->setRange(0,1.0);
  R->setBins(25);

  TString fullname = Form("%s_%s",name.Data(),cats.at(catIndex).Data());

  TString selectionString = "";

  fitInfo *info = & (per_cat_fit_ranges.at(catIndex));

  switch(reg) {
  case kSignal:
    fullname+="_Signal";
    selectionString = Form("mgg > %f && mgg < %f",info->signal_min,info->signal_max);
    break;
  case kSideband:
    fullname+="_Sideband";
    selectionString = Form("(mgg > %f && mgg < %f) || (mgg > %f && mgg < %f)",
			   info->sideband_low_min,  info->sideband_low_max,
			   info->sideband_high_min, info->sideband_high_max);
    break;
  case kAll:
    fullname+="_All";
    selectionString = "1";
    break;
  default:
    throw new std::runtime_error( Form("invalid region %d",reg) );
  }


  TH2F* hist = ((RooDataSet*)ws->data(Form("%s_%s",name.Data(),cats.at(catIndex).Data()))->reduce(selectionString))->createHistogram(*MR,*R,"",fullname);
  hist->SetName(fullname);
  hist->Scale(weight);

  RooDataHist roohist(Form("%s_hist",fullname.Data()),"",RooArgSet(*MR,*R),*((RooDataSet*)ws->data(Form("%s_%s",name.Data(),cats.at(catIndex).Data()))->reduce(selectionString)),
		   weight);
  RooHistPdf roopdf(Form("%s_pdf",fullname.Data()),"",RooArgSet(*MR,*R),roohist);

  std::cout << roohist.GetName() << std::endl;

  ws->import(roohist);
  ws->import(roopdf);


  return hist;

}

void Fitter::saveAll() {
  outputFile->cd();
  ws->Write();
  
  for(auto& it: save_histograms) {
    it.second->Write();
    getSignalRegionHistogram(*it.second,it.first+"_SigRegions")->Write();
  }
}

void Fitter::defineCats() {
  cats.push_back("SingleMu");
  cats.push_back("SingleEle");
  cats.push_back("Hadronic_EBEB");

  per_cat_fit_ranges.push_back({117.,121.,123.5,127.5,129.,133.});
  per_cat_fit_ranges.push_back({117.,121.,123.5,127.5,129.,133.});
  per_cat_fit_ranges.push_back({117.,121.,123.5,127.5,129.,133.});

  sideband_integrals.resize(cats.size());
  signal_integrals.resize(cats.size());
  sideband_to_signal_scale_factors.resize(cats.size());
  setupCategory("data");
  for(auto& it: sm_higgs_names) {
    setupCategory(it);
  }
  for(auto& it: sms_names) {
    setupCategory(it);
  }
}

void Fitter::setupCategory(TString name) {
  {
    RooDataSet *data = (RooDataSet*)ws->data(name);
    data = (RooDataSet*)data->reduce("mu1_pt > 20.");
    ws->import(*data,RooFit::Rename(Form("%s_SingleMu",name.Data())) );
  }
  
  {
    RooDataSet *data = (RooDataSet*)ws->data(name);
    data = (RooDataSet*)data->reduce("mu1_pt <= 20. && ele1_pt > 20.");
    ws->import(*data,RooFit::Rename(Form("%s_SingleEle",name.Data())) );
  }

  {
    RooDataSet *data = (RooDataSet*)ws->data(name);
    data = (RooDataSet*)data->reduce("mu1_pt <= 20. && ele1_pt <=20.");
    ws->import(*data,RooFit::Rename(Form("%s_Hadronic_EBEB",name.Data())) );
  }

}
