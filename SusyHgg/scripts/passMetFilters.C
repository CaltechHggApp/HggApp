#include "TTree.h"

#include "TFile.h"
#include "TString.h"

#include <fstream>
#include <unordered_set>
#include <vector>
#include "stdint.h"

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

  Long64_t iEntry=-1;
  while(analysisTree->GetEntry(++iEntry)) {
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
    "HBHENoiseFilterFlag",
    "isNotDeadEcalCluster",
    "trackerFailureFilterFlag",
    "CSCHaloFilterFlag",
    "ECALTPFilterFlag"
  };
  uint32_t *flags = new uint32_t[metFlags.size()];

  for(int i=0;i<metFlags.size();i++) {
    reducedTree->SetBranchStatus(metFlags.at(i),1);
    reducedTree->SetBranchAddress(metFlags.at(i),&flags[i]);
  }

  iEntry=-1;
  while(reducedTree->GetEntry(++iEntry)) {
    if(selected.count(evtInfo(rRun,rLumi,rEvt))==0) continue;
    bool pass=1;
    for(int i=0;i<metFlags.size();i++) pass = (pass && flags[i]);
    if(pass) passing.insert(evtInfo(rRun,rLumi,rEvt));
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
