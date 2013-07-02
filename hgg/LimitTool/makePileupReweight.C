#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TString.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
//using namespace std;

void makePileupReweight(TChain *chain, TH1D *dataHist,TString outputName){
  std::cout << "makePileupReweight" << std::endl;
  TH1D* mcPU = new TH1D("mcPU","",dataHist->GetNbinsX(),0,dataHist->GetNbinsX());
  std::cout << "project" << std::endl;
  chain->Project("mcPU","nPU[1]");
  std::cout << "project" << std::endl;
  mcPU->Scale(1./mcPU->Integral());
  
  dataHist->Scale(1./dataHist->Integral());
  std::cout << "clone" << std::endl;
  TH1D* pu = (TH1D*)dataHist->Clone("pileupReWeight");
  pu->Divide(mcPU);
  std::cout << "clone" << std::endl;
  
  std::cout << "write"<< std::endl;
  TFile *f = new TFile(outputName,"RECREATE");
  pu->Write();
  f->Close();
  std::cout << "write"<< std::endl;
  delete mcPU;
  delete pu;
}

void makeAllPU(string mcFileList,TString dataFileName,TString outputFile){
  TFile * dataFile = new TFile(dataFileName);
  TH1D* dataHist = (TH1D*)dataFile->Get("pileup");
  dataHist->Scale(1./dataHist->Integral());

  string mcFile;
  ifstream *s = new ifstream(mcFileList.c_str());
  TChain* chain = new TChain("ntp1");
  while(s->good()){
    getline(*s,mcFile);
    std::cout << mcFile << std::endl;
    if(mcFile.compare("")==0) continue;
    chain->AddFile(mcFile.c_str());
  }
  makePileupReweight(chain,dataHist,outputFile);
  dataFile->Close();
  s->close();
  delete chain;

}

void make1PU(string mcFileList,TString dataFileName,TString outputFile){
  TFile * dataFile = new TFile(dataFileName);
  TH1D* dataHist = (TH1D*)dataFile->Get("pileup");
  dataHist->Scale(1./dataHist->Integral());

  string mcFile;
  ifstream *s = new ifstream(mcFileList.c_str());
  TChain *chain = new TChain("HggReduce");
  while(s->good()){
    getline(*s,mcFile);
    chain->Add(mcFile.c_str());
  }
  makePileupReweight(chain,dataHist,outputFile);
  dataFile->Close();
  s->close();

}

