#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <cstdint>
#include "string.h"

#include "ArgParser.hh"
#include "assert.h"

using namespace std;

struct evtInfo { 
  uint32_t run,lumi,event;  
  evtInfo(uint32_t r, uint32_t l, uint32_t e){run=r; lumi=l; event=e;}
  evtInfo(int32_t r, int32_t l, int32_t e){run=static_cast<uint32_t>(r); lumi=static_cast<uint32_t>(l); event=static_cast<uint32_t>(e);}
 bool operator== (const evtInfo& e) const{
    return ((e.run==run) && (e.lumi==lumi) && (e.event==event));
  }

};

struct eiHash{
  size_t operator() (const evtInfo& ei) const {
    return ei.run+size_t(ei.lumi)<<18+size_t(ei.event)<<31;
  }
};

void passMetFilters(TTree* reducedTree,TTree* analysisTree,const char* outputFile) {

  std::unordered_set<evtInfo,eiHash> selected;
  std::unordered_set<evtInfo,eiHash> passing;

  int aRun,aLumi,aEvt;
  

  analysisTree->SetBranchAddress("run",&aRun);
  analysisTree->SetBranchAddress("lumi",&aLumi);  
  analysisTree->SetBranchAddress("evt",&aEvt);
  
  std::cout << "First Pass Analysis Tree" << std::endl;
  Long64_t iEntry=-1;
  while(analysisTree->GetEntry(++iEntry)) {
    if(iEntry%5000==0) std::cout << "Processing Entry " << iEntry << std::endl;
    selected.insert(evtInfo(aRun,aLumi,aEvt));
  }

  int rRun,rLumi,rEvt;
  
  reducedTree->SetBranchStatus("*",0);
  reducedTree->SetBranchStatus("runNumber",1);
  reducedTree->SetBranchStatus("lumiBlock",1);
  reducedTree->SetBranchStatus("evtNumber",1);

  reducedTree->SetBranchAddress("runNumber",&rRun);
  reducedTree->SetBranchAddress("lumiBlock",&rLumi);
  reducedTree->SetBranchAddress("evtNumber",&rEvt);
  
  std::vector<TString> metFlags = {
    "eeBadScFilterFlag",
    "hcalLaserEventFilterFlag",
    "HBHENoiseFilterResultFlag",
    //"isNotDeadEcalCluster",
    "trackerFailureFilterFlag",
    "CSCHaloFilterFlag",
    "ECALTPFilterFlag"
  };
  uint8_t *flags = new uint8_t[metFlags.size()];

  for(int i=0;i<metFlags.size();i++) {
    reducedTree->SetBranchStatus(metFlags.at(i),1);
    reducedTree->SetBranchAddress(metFlags.at(i),&flags[i]);
  }

  iEntry=-1;
  uint32_t nFailingNoise=0;
  while(reducedTree->GetEntry(++iEntry)) {
    if(iEntry%5000==0) std::cout << "Processing Entry " << iEntry << "   " << passing.size() << " events selected   " << nFailingNoise << " failing noise filter" << std::endl;
    if(selected.count(evtInfo(rRun,rLumi,rEvt))==0) continue;
    bool pass=1;
    for(int i=0;i<metFlags.size();i++) {
      pass = (pass && flags[i]);
      if(flags[i]!=0 && flags[i]!=1) {
	std::cout << "ERROR:  " << metFlags.at(i) << ":  " << uint16_t(flags[i]) << std::endl;
	assert(false);
      }
      if(rRun==191833 && rEvt==86217002) std::cout << "\t" << metFlags.at(i) << " :  " << uint16_t(flags[i]) << "      " << pass << std::endl;
    }
    if(pass) passing.insert(evtInfo(rRun,rLumi,rEvt));
    else nFailingNoise++;
  }

  TFile outF(outputFile,"RECREATE");
  TTree metNoiseTree("SusyHggNoiseTree","");
  bool passFilters;
  metNoiseTree.Branch("passNoise",&passFilters,"passNoise/O");
  iEntry=-1;
  while(analysisTree->GetEntry(++iEntry)) {
    passFilters = (passing.count(evtInfo(aRun,aLumi,aEvt))==1);
    metNoiseTree.Fill();
  }
  metNoiseTree.Write();
  outF.Close();
}


int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("ReducedTree",ArgParser::required,"List of Reduced nTuples");
  a.addArgument("AnalysisTree",ArgParser::required,"path to analysis nTuple");
  a.addArgument("OutputFile",ArgParser::required,"output file");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  TChain *theChain = new TChain("HggReduce");
  std::cout << "chaining" << std::endl;
  char Buffer[2000];
  TString RootFileName;
  char MyRootFile[2000];  
  ifstream *inputFile = new ifstream(a.getArgument("ReducedTree").c_str());
  if(!inputFile){
    std::cout << "Invalid input file" <<std::endl;
    return 1;
  }
  // get the tree with the conditions from the first file
  //  TTree *treeCond = new TTree();
  //  int nfiles=1;
  
  bool useXRD = a.longFlagPres("XRootD");
  TString xrdSite="";
  if(useXRD) xrdSite = a.getLongFlag("XRootD");
  
  char tmpFileName[2000];
  
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,2000);
    RootFileName = TString(Buffer).Strip();
    
    if(RootFileName.First('#') != -1 || RootFileName.Length()==0) continue;
    
    if(useXRD) {
      RootFileName.Replace(0,RootFileName.Index("/store"),
			   Form("root://xrootd.unl.edu//store/test/xrootd/%s/",xrdSite.Data()));
    }
    theChain->Add(RootFileName);
    
    std::cout << "chaining " << RootFileName << std::endl;
  } 
  inputFile->close();
  
  TFile * analysisTreeFile = new TFile(a.getArgument("AnalysisTree").c_str());
  TTree * analysisTree     = (TTree*)analysisTreeFile->Get("SusyHggTree");

  

  passMetFilters(theChain,analysisTree,a.getArgument("OutputFile").c_str());
  
}
