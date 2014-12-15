#include "TTree.h"
#include "TString.h"
#include "TFile.h"

#include <iostream>

float getSherpaWeight(int nJets) {
  if(nJets<3) return 1.;
  if(nJets==3) return 0.5;
  if(nJets==4) return 0.2;
  return 0.1;
}

void weightBkgMC(TString inputFileName,TString outputFileName,float weight,TString treeName="SusyHggTree") {
  
  float inw,outw;
  
  bool isSherpa = (inputFileName.Index("DiPhotonJetsBox_M60_8TeV-sherpa") != -1);


  int Nj;

  TFile inputFile(inputFileName);
  TTree* inputTree = (TTree*)inputFile.Get(treeName);

  TFile outputFile(outputFileName,"RECREATE");
  TTree* outputTree = inputTree->CloneTree(0);

  inputTree->SetBranchAddress("pileupWeight",&inw);
  inputTree->SetBranchAddress("Njets",&Nj);
  outputTree->SetBranchAddress("pileupWeight",&outw);

  Long64_t iEntry=-1;

  if(isSherpa) {
    double wSum=0;
    double entries=0;
    while(inputTree->GetEntry(++iEntry)) {
      wSum+=getSherpaWeight(Nj);
      entries++;
    }
    weight = weight*entries/wSum;

    std::cout << "sherpa: sum of weights: " << wSum << "   total entries: " << entries << std::endl;
  }
  iEntry=-1;
  while(inputTree->GetEntry(++iEntry)) {
    outw = inw*weight;
    if(isSherpa) outw*=getSherpaWeight(Nj);
    outputTree->Fill();
  }

  outputFile.cd();
  outputTree->Write();
  outputFile.Close();
  inputFile.Close();

}

