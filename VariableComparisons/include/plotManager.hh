#ifndef plotManager_hh
#define plotManager_hh

#include "TString.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeFormula.h"
#include "TFile.h"

#include <array>
#include <map>
#include <vector>
#include <iostream>

#include "include/varCorrector.hh"

class plotManager{
public:
  plotManager(TString inTag="");
  void addCategory(TString name,TString cut);
  void addVariable(TString name,TString var,int bins, float low, float high,bool correct=false);
  void setTargetPU(TH1D* hist) { target_pu=hist;}
  void setMCPU(TH1D* hist) { mc_pu=hist;}
  virtual void processChain(TChain *fChain,float weight);
  std::array<TH1F*,3> getHistogram(TString cat, TString var);
  void addVetos(std::vector<TString>* v){vetos = v;}
  void saveAll(TFile *f);
  void setUse4Cat(bool b = true){use4Cat=b;}


  void buildHistograms();
protected:
  std::vector<TString>* vetos;
  std::vector<TTreeFormula*> vetoFormulas;
  bool freeze;
  TString histNameTag;
  bool use4Cat;
  varCorrector corrector;
  varCorrector4cat corrector4cat;

  std::vector<TString> catNames;
  std::vector<TString> catCuts;
  TTreeFormula* etaVal;
  TTreeFormula* isRealPho;
  TTreeFormula* isRealEle;
  TTreeFormula* isTrigger;
  TTreeFormula* pu;


  std::vector<TString> varNames;
  std::vector<TString> variables;
  std::vector<int> xBins;
  std::vector<float> xLow;
  std::vector<float> xHigh;
  std::vector<bool> useCorrection;

  std::vector< std::vector< TH1F* > > fakeHistograms; // outer vector: category % inner vector: variables
  std::vector< std::vector< TH1F* > > realPhoHistograms; // outer vector: category % inner vector: variables
  std::vector< std::vector< TH1F* > > realEleHistograms; // outer vector: category % inner vector: variables
  std::vector< TTreeFormula* > catFormulas;
  std::vector< TTreeFormula* > varFormulas;  
  void buildFormulas(TChain *fChain);
  void updateFormulas();
  void destroyFormulas();
  void processEntry(float weight);

  float computePUWeight();

  TH1D* target_pu=0;
  TH1D* mc_pu=0;
};



#endif
