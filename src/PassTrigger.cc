#include "PassTrigger.hh"

#include <fstream>
#include <iostream>

#include <stdio.h>

PassTrigger::PassTrigger(TTree* vecbostree,TTree* analysistree) : Vecbos(vecbostree) {
  analysisTree=analysistree;

  analysisTree->SetBranchAddress(runLabel.c_str(),&run);
  analysisTree->SetBranchAddress(lumiLabel.c_str(),&lumi);  
  analysisTree->SetBranchAddress(evtLabel.c_str(),&evt);
}

void PassTrigger::Loop(string outputFileName) {
  if(outputStep){
    LoopOnOutput(outputFileName);
  }else{
    LoopOnVecbos(outputFileName);
  }
}


void PassTrigger::LoopOnOutput(string outputFileName) {
  //read the text file
  char tmpFileName[L_tmpnam];
  tmpnam(tmpFileName);
  TFile tmpFile(tmpFileName,"RECREATE");
  TTree tmp("tmp","");
  tmp.ReadFile(inputText.c_str());
  int r,l,e;
  tmp.SetBranchAddress("run",&r);
  tmp.SetBranchAddress("lumi",&l);
  tmp.SetBranchAddress("evt",&e);
  Long64_t iEntry=-1;
  while(tmp.GetEntry(++iEntry)) {
    if(e==480229083 || e==50763554) std::cout << r << " " << l << " " << e << std::endl;
    passTriggerEvents.insert(evtInfo(static_cast<uint32_t>(r),static_cast<uint32_t>(l),static_cast<uint32_t>(e)));
  }
  tmpFile.Close();
  
  //make the output tree
  std::cout << "Found " << passTriggerEvents.size() << " events passing the trigger" << std::endl;
  TFile outputRootFile(outputFileName.c_str(),"RECREATE");
  iEntry=-1;
  TTree outputTree("SusyHggTriggerTree","");
  bool passHLT;
  int nPass=0;
  outputTree.Branch("passTrigger",&passHLT,"passTrigger/O");
  while(analysisTree->GetEntry(++iEntry)) {
    if(iEntry%5000==0) cout << "Processing Entry " << iEntry << "\t\t" << nPass << " events passed trigger" << std::endl;
    passHLT = (passTriggerEvents.count(evtInfo(static_cast<uint32_t>(run),static_cast<uint32_t>(lumi),static_cast<uint32_t>(evt)))==1);
    if(evt==480229083 || evt==50763554) std::cout << evt << "  " << passHLT << std::endl;
    nPass+=passHLT;
    outputTree.Fill();
  }
  outputTree.Write();
  outputRootFile.Close();
}

void PassTrigger::LoopOnVecbos(string outputFileName) {
  Long64_t iEntry =-1;
    cout << "Processing Analysis Tree" << std::endl;
    while(analysisTree->GetEntry(++iEntry)) {
      if(iEntry%5000==0) cout << "Processing Entry " << iEntry << std::endl;
      selectedEvents.insert(evtInfo(static_cast<uint32_t>(run),static_cast<uint32_t>(lumi),static_cast<uint32_t>(evt)));
    }
    cout << "Done\n\tSelected " << selectedEvents.size() << " events" << std::endl;
    
    cout << "Processing vecbos tree" << std::endl;
    ofstream outFile(outputFileName);
    outFile << "run/I:lumi/I:evt/I\n";
    
    iEntry=-1;
    
    fChain->SetBranchStatus("*",0);
    fChain->SetBranchStatus("runNumber",1);
    fChain->SetBranchStatus("eventNumber",1);
    fChain->SetBranchStatus("lumiBlock",1);
    fChain->SetBranchStatus("*Trg",1);
    fChain->SetBranchStatus("*HLT",1);
    setRequiredTriggers(triggers);
    int nFound=0;
    while(fChain->GetEntry(++iEntry)) {
      if(iEntry %5000==0) cout << "Processing Entry " << iEntry << std::endl;
      if(selectedEvents.count(evtInfo(static_cast<uint32_t>(runNumber),static_cast<uint32_t>(lumiBlock),static_cast<uint32_t>(eventNumber)))==0) continue;
      nFound++;
      reloadTriggerMask(true);
      if(hasPassedHLT()) {
	outFile << runNumber << " " << lumiBlock << " " << eventNumber << std::endl;
	passTriggerEvents.insert(evtInfo(static_cast<uint32_t>(runNumber),static_cast<uint32_t>(lumiBlock),static_cast<uint32_t>(eventNumber)));
      }
    }
    outFile.close();
}
