#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TTreeFormula.h"
#include "TLatex.h"
#include "TList.h"

#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <array>
#include <assert.h>
#include <string>
#include <tuple>

#include "../include/ArgParser.hh"
#include "include/HistogramStack.hh"


#include <cstdlib>
#include <algorithm>
#include <iomanip>


//local includes
#include "include/weightManager.hh"
#include "include/varCorrector.hh"
#include "include/plotManager2D.hh"
#include "src/getCategories.C"


//TCanvas *makeCanvas(std::array<TH1F*,3> data,std::array<TH1F*,3> mc, TString xName,TString label="");

void writeRootFile(TString fileName,bool useSherpa,bool test);
void CalcFromRootFile(TString fileName,bool useSherpa,bool test);

int main(int argc,char** argv){
  ArgParser a(argc,argv);


  a.addArgument("mode",ArgParser::required,"specify the operation mode (save,calculate)");
  a.addArgument("rootFileName",ArgParser::required,"specify the name of the root file");
  a.addLongOption("sherpa",ArgParser::noArg,"use SHERPA diphoton jets sample");
  a.addLongOption("test",ArgParser::noArg,"Only process 1 MC file (useful for testing)");
  a.addShortOption('h',ArgParser::noArg,"print help");
  std::string ret;
  if(a.process(ret)!=0) {
    std::cout << "\nInvalid Option " << ret << std::endl << std::endl;;
    a.printOptions(argv[0]);
    return 0;
  }

  if(a.shortFlagPres('h')) {
    a.printOptions(argv[0]);
    return 0;
  }
  
  enum OpMode{kNone,kSave,kCalc};
  OpMode mode = kNone;

  if(a.getArgument("mode").compare("save")==0) mode=kSave;
  if(a.getArgument("mode").compare("calculate")==0) mode=kCalc;
  
  if(mode==kNone) {
    std::cout << "\nInvalid Mode Speficied:  " << a.getArgument("mode") << std::endl << std::endl;
    return -1;
  }
  TString fileName = a.getArgument("rootFileName");

  bool useSherpa = a.longFlagPres("sherpa");
  bool isTest    = a.longFlagPres("test");

  switch(mode) {
  case kSave:
    writeRootFile(fileName,useSherpa,isTest);
    break;
  case kCalc:
    CalcFromRootFile(fileName,useSherpa,isTest);
    break;
  }


  return 0;
}

void writeRootFile(TString fileName,bool useSherpa,bool test) {
  if(useSherpa) std::cout << "SHERPA!!" << std::endl;

  float lumi=6.3;
  weightManager weights;

  plotManager2D mcPlotter("mc");
  plotManager2D dataPlotter("data");
  mcPlotter.setUse4Cat();
  dataPlotter.setUse4Cat();

  std::vector<TString> samples = {
    "GJets_Pt20to40.root",
    "GJets_Pt40.root",
    "QCD_Pt30to40.root",
    "QCD_Pt40.root",
    "DYJetsToLL_M-50.root" };


  if(useSherpa) {
    samples.push_back("DiPhotonJets_sherpa.root");
  }
  else {
    samples.push_back("DiPhotonBox_Pt10to25.root");
    samples.push_back("DiPhotonBox_Pt25to250.root");
    samples.push_back("DiPhotonBox_Pt250.root");
    samples.push_back("DiPhotonJets.root");    
  }

  int Nsamples = samples.size();

  std::vector<TString> data = {
    "DoublePhoton_Run2012B_13Jul2012.root"
  };

  int Ndata = data.size();


  std::vector<TString> globalVetos;
  globalVetos.push_back("etaSC > -1.775 && etaSC < -1.764 && phi > 1.2 && phi < 1.6");
  globalVetos.push_back("etaSC > 1.591 && etaSC < 1.595 && phi > -2.06 && phi < -2.045");
  //globalVetos.push_back("etaSC > -0.9 && etaSC < -0.7 && phi > 2.9 && phi < 3.1");
  globalVetos.push_back("etaSC > 1.74 && etaSC < 1.76 && phi > 2.1 && phi < 2.15");
  globalVetos.push_back("etaSC > 1.564 && etaSC < 1.565 && phi > 0.528 && phi < 0.532");

  //output->Draw("phi:etaSC","!(etaSC > -1.8 && etaSC < -1.76 && phi > 1.2 && phi < 1.5) && !(etaSC>1.591 && etaSC<1.595 && phi > -2.06 && phi < -2.045) && !(etaSC > 1.74 && etaSC < 1.76 && phi > 2.1 && phi < 2.15) && !(etaSC > 1.564 && etaSC < 1.565 && phi > 0.528 && phi < 0.532) && !(etaSC > -1.6888 && etaSC < -1.6868 && phi > 1.8 && phi < 1.802) && sieie < 1e-5 && abs(etaSC) > 1.557 && abs(etaSC)<2

  auto catInfo = getCategories();
  auto cats = catInfo.first;
  auto catNames = catInfo.second;



  assert(cats.size() == catNames.size());

  int Ncat = cats.size();

  if(test) {
    Ncat = Ndata = Nsamples = 1;
  }

  for(int i=0;i<Ncat;i++){
    mcPlotter.addCategory(catNames[i],cats[i]);
    dataPlotter.addCategory(catNames[i],cats[i]);
  }

  mcPlotter.addVetos(&globalVetos);
  dataPlotter.addVetos(&globalVetos);

  auto vars = getVars();
  for(auto it : vars) {
    mcPlotter.addVariable  (translateVar(std::get<0>(it)),std::get<0>(it),std::get<1>(it),std::get<2>(it),std::get<3>(it));
    dataPlotter.addVariable  (translateVar(std::get<0>(it)),std::get<0>(it),std::get<1>(it),std::get<2>(it),std::get<3>(it));
  }

  int Nvar = vars.size();

  for(int iSample=0;iSample<Nsamples;iSample++){
    TFile *f = new TFile("output/"+samples[iSample]);
    TChain *fChain = (TChain*)f->Get("output");
    TString sampleName = samples[iSample];
    sampleName.Remove(sampleName.Last('.'));
    std::cout << sampleName << std::endl;
    mcPlotter.processChain(fChain,weights.getWeight(sampleName,lumi));
  }

  for(int iData=0;iData<Ndata;iData++){
    TFile *f = new TFile("output/"+data[iData]);
    TChain *fChain = (TChain*)f->Get("output");
    dataPlotter.processChain(fChain,1);
  }


  TFile outputFile(fileName,"RECREATE");
  mcPlotter.saveAll2D(&outputFile);
  dataPlotter.saveAll2D(&outputFile);

} //void writeRootFile

struct plotRef{
  TString var1,var2,cat;
  float mc_corr,data_corr;
  float diff;
  void calcDiff() { diff = fabs(mc_corr-data_corr); }
};


void getCorrFromRef(plotRef& ref,TFile *file) {
  std::array<TString,3> histSuffix = {"realPho","realEle","fake"};
  TH2F* mc =   (TH2F*)file->Get( ref.var1+"_"+ref.var2+"_"+ref.cat+"_mc_"+histSuffix[0] );
  if(mc==0) return;
  TH2F* data = (TH2F*)file->Get( ref.var1+"_"+ref.var2+"_"+ref.cat+"_data_"+histSuffix[0] );

  for(int i=1;i<3;i++) {
    mc->Add(  (TH2F*)file->Get( ref.var1+"_"+ref.var2+"_"+ref.cat+"_mc_"+histSuffix[i] )  );
    data->Add(  (TH2F*)file->Get( ref.var1+"_"+ref.var2+"_"+ref.cat+"_data_"+histSuffix[i] )  );
  }

  ref.mc_corr = mc->GetCorrelationFactor();
  ref.data_corr = data->GetCorrelationFactor();
  ref.calcDiff();
}

bool ComparePlotRef(const plotRef& a, const plotRef& b) {
  return a.diff > b.diff;
  //return a.data_corr < b.data_corr;
}

void CalcFromRootFile(TString fileName,bool useSherpa,bool isTest) {
  auto catInfo = getCategories();
  auto vars     = getVars();


  std::vector<TString> catNames = catInfo.second;
  int Ncat = catNames.size();
  int Nvar = vars.size();

  TFile *file = new TFile(fileName);

  std::vector<plotRef> correlationList;

  std::array<TString,3> histSuffix = {"realPho","realEle","fake"};
  for(auto catIt : catNames) {
    for(auto var1It : vars) {
      for(auto var2It: vars) {
	plotRef r = {std::get<0>(var1It),std::get<0>(var2It),catIt,-2,2};
	getCorrFromRef(r,file);
	if(r.mc_corr<-1) continue;
	correlationList.push_back(r);
      }
    }
  }
  std::sort(correlationList.begin(),correlationList.end(),ComparePlotRef);
  
  //std::cout << "var1    var2   cat     corr_mc     corr_data      diff" << std::endl;
  std::cout << std::setw(50) << "var1";
  std::cout << std::setw(50) << "var2";
  std::cout << std::setw(80) << "cat";
  std::cout << std::setw(25) << "Corr MC";
  std::cout << std::setw(25) << "Corr Data";
  std::cout << std::endl;
  int i=0;
  for(auto r: correlationList) {
    if(++i >=60) break;
  std::cout << std::setw(50) << r.var1;
  std::cout << std::setw(50) << r.var2;
  std::cout << std::setw(80) << r.cat;
  std::cout << std::setw(25) << r.mc_corr;
  std::cout << std::setw(25) << r.data_corr;
  std::cout << std::endl;
  }
}
