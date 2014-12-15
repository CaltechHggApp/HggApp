#include "TChain.h"
#include "TString.h"
#include "TTreeFormula.h"
#include "TFile.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooFitResult.h"

#include <map>

TTreeFormula* getBaseSelection(TChain* tree, bool isDY);

TTreeFormula* getEta(TChain* tree, bool isDY);
TTreeFormula* getR9(TChain* tree, bool isDY);
TTreeFormula* getN(TChain* tree, bool isDY);

TTreeFormula* getDeltaENoReg(TChain* tree, bool isDY);
TTreeFormula* getDeltaE(TChain* tree, bool isDY);
TTreeFormula* getEP(TChain* tree, bool isDY);

//-1=veto, 0-EBLow HighR9, 1-EBHigh High R9, ..., 4- EBLow LowR9,...
int getEtaR9Region(float eta, float r9); 
float getDeltaEScale(int region);

RooFitResult* fitCB(RooDataSet* data, RooRealVar* var);

void MeasureSmearing(TChain* tree, TString outputTag,bool isDY=false,bool doScaling=false) {

  std::map<TString,TTreeFormula*> formulas;

  formulas["base"] = getBaseSelection(tree,isDY);
  formulas["Eta"] = getEta(tree,isDY);
  formulas["R9"] = getR9(tree,isDY);
  formulas["N"] = getN(tree,isDY);
  formulas["DeltaENoReg"] = getDeltaENoReg(tree,isDY);
  formulas["DeltaE"] = getDeltaE(tree,isDY);
  formulas["EP"] = getEP(tree,isDY);

  RooRealVar denr("DeltaENoReg","",-0.5,1.5);
  RooRealVar de("DeltaE","",-0.5,1.5);
  RooRealVar ep("EP","",0,3);

  std::vector<RooDataSet*> datasets;
  for(int i=0;i<8;i++) datasets.push_back( new RooDataSet( Form("data_%d",i),"",RooArgSet(denr,de,ep)) );
  
  tree->SetBranchStatus("*",0);
  if(isDY) {
    tree->SetBranchStatus("nEle",1);
    tree->SetBranchStatus("Electrons*",1);
  }else{
    tree->SetBranchStatus("nPho",1);
    tree->SetBranchStatus("Photons*",1);
  }    

  int iTree=-1;
  Long64_t iEntry=-1;
  while(tree->GetEntry(++iEntry)) {
    if(iEntry%1000==0) std::cout << "Processing Entry " << iEntry << std::endl;
    //update formulas if necessary
    if(tree->GetTreeNumber() != iTree) {
      //update all the TTreeFormulas
      for(std::map<TString,TTreeFormula*>::iterator fIt = formulas.begin();
	  fIt!= formulas.end(); fIt++) {
	if(fIt->second) fIt->second->UpdateFormulaLeaves();       
      }
      //change the tree number
      iTree = tree->GetTreeNumber();
    }
    int N = formulas["N"]->EvalInstance();
    for(int iPho=0;iPho<N;iPho++) {

      if(formulas["base"]->EvalInstance(iPho) == false) continue;
      int cat = getEtaR9Region(formulas["Eta"]->EvalInstance(iPho), formulas["R9"]->EvalInstance(iPho));
      if(cat<0) continue;
      
      float dENR = formulas["DeltaENoReg"]->EvalInstance(iPho);
      float dE = formulas["DeltaE"]->EvalInstance(iPho);
      float EP = 0;
      if(isDY) EP = formulas["EP"]->EvalInstance(iPho);

      denr.setVal(dENR);
      de.setVal(dE);
      ep.setVal(EP);

      datasets.at(cat)->add(RooArgSet(denr,de,ep));
    }

  }

  std::vector<RooFitResult*> fitDENR,fitDE, fitEP;
  for(int i=0;i<8;i++) {
    fitDENR.push_back( (RooFitResult*)fitCB( datasets.at(i), &denr )->Clone(Form("fitres_denr_%d",i)));
    fitDE.push_back((RooFitResult*)fitCB( datasets.at(i), &de )->Clone(Form("fitres_de_%d",i)));
    if(isDY) fitEP.push_back((RooFitResult*)fitCB( datasets.at(i), &ep )->Clone(Form("fitres_ep_%d",i)));
  }

  for(int i=0;i<8;i++) {
    std::cout << "category " << i << std::endl;
    std::cout << "DE No Reg" << std::endl;
    fitDENR.at(i)->Print();
    fitDENR.at(i)->Print("V");
    std::cout << "DE" << std::endl;
    fitDE.at(i)->Print();
    fitDE.at(i)->Print("V");
    std::cout << "EP" << std::endl;
    fitEP.at(i)->Print();
    fitEP.at(i)->Print("V");
  }

  //delete formulas
  for(std::map<TString,TTreeFormula*>::iterator fIt = formulas.begin();
      fIt!= formulas.end(); fIt++) {
    delete fIt->second;
  }
}


TTreeFormula* getBaseSelection(TChain* tree, bool isDY) {
  TString sel = "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25";
  if(isDY) sel.ReplaceAll("Photons","Electrons");
  return new TTreeFormula("baseSelection",sel,(TTree*)tree);
}

TTreeFormula* getEta(TChain* tree, bool isDY) {
  TString sel = "Photons.SC.eta";
  if(isDY) sel.ReplaceAll("Photons","Electrons");
  return new TTreeFormula("eta",sel,(TTree*)tree);
}

TTreeFormula* getR9(TChain* tree, bool isDY) {
  TString sel = "Photons.SC.r9";
  if(isDY) sel.ReplaceAll("Photons","Electrons");
  return new TTreeFormula("r9",sel,(TTree*)tree);
}

TTreeFormula* getN(TChain* tree, bool isDY) {
  TString sel = "nPho";
  if(isDY) sel.ReplaceAll("Pho","Ele");
  return new TTreeFormula("N",sel,(TTree*)tree);
}

TTreeFormula* getDeltaENoReg(TChain* tree, bool isDY) {
  TString sel = "Photons.energy/Photons.genMatch.energy-1";
  if(isDY) sel.ReplaceAll("Photons","Electrons");
  return new TTreeFormula("DeltaE",sel,(TTree*)tree);  
}

TTreeFormula* getDeltaE(TChain* tree, bool isDY) {
  TString sel = "Photons.correctedEnergy/Photons.genMatch.energy-1";
  if(isDY) sel.ReplaceAll("Photons","Electrons");
  return new TTreeFormula("DeltaE",sel,(TTree*)tree);  
}

TTreeFormula* getEP(TChain* tree, bool isDY) {
  if(!isDY) return 0;
  TString sel = "Electrons.EOverP";
  return new TTreeFormula("EP",sel,(TTree*)tree);  
}



int getEtaR9Region(float eta,float r9) {
  if( fabs(eta)>1.442 && fabs(eta)<1.56) return -1;
  if(fabs(eta)>2.5) return -1;

  int reg=0;
  if(r9<0.94) reg+=4;

  if(fabs(eta)>=1.) reg++;
  if(fabs(eta)>=1.56) reg++;
  if(fabs(eta)>=2.) reg++;

  return reg;
}

float getDeltaEScale(int region) {

  if(region==2) return 1.0276; //EElow high r9
  if(region==3) return 1.0216; //EEhigh high r9

  if(region==6) return 1.0287; //EElow low r9
  if(region==7) return 1.0351;

  return 1;
}

RooFitResult* fitCB(RooDataSet* data, RooRealVar* var) {
  RooRealVar peak("peak","",-1,3);
  RooRealVar sig("sig","",0.01,4);
  RooRealVar alpha("alpha","",0.001,10);
  RooRealVar n("n","",0.001,5);

  RooCBShape cb("cb","",*var,peak,sig,alpha,n);

  cb.fitTo(*data,RooFit::Save(kFALSE),RooFit::Strategy(0));
  RooFitResult* res = cb.fitTo(*data,RooFit::Save(kFALSE),RooFit::Strategy(2));

  return res;
}
