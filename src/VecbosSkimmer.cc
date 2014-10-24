#include "VecbosSkimmer.hh"
#include "TTree.h"
#include <iostream>

void VecbosSkimmer::Skim(const char* outputFileName) {
  int64_t iEntry=-1;

  TFile * outputFile = new TFile(outputFileName,"RECREATE");
  TTree * outputTree = fChain->CloneTree(0);

  uint32_t nSel=0;

  while(fChain->GetEntry(++iEntry)) {
    if(nSel >= eventsToSkim.size()) break;
    if(iEntry % 500 == 0) std::cout << "Processing " << iEntry << "   " 
				     << nSel << " / " << eventsToSkim.size() << "  Selected" << std::endl;
    if( eventsToSkim.count(evtInfo(runNumber,lumiBlock,eventNumber)) == 1) {
      outputTree->Fill();
      nSel++;
    } 
  }

  outputTree->Write();
  outputFile->Close();
}
