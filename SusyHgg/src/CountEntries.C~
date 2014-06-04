#include <iostream>
#include <string>
#include <fstream>

#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TTreeFormula.h"

#include "../include/ArgParser.hh"

int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("InputList","file listing the input .root files");
  a.addArgument("OutputFileName","Name of the Output File");

  std::string ret;
  if(a.process(ret) !=0){
    std::cout << "Invalid Options:  " << ret << std::endl;
    a.printOptions(argv[0]);
    return 0;
  }

  std::string name = "ntp1";
  TChain c(name.c_str(),name.c_str());

  std::string input = a.getArgument("InputList");
  std::string output = a.getArgument("OutputFileName");
  fstream listfile(input.c_str(),std::fstream::in);
  char buffer[2000];

  while(listfile.getline(buffer,2000)) {
    c.Add(buffer);
  }

  TH2F N("Ntotal","",20,12.5,512.5,20,12.5,512.5);
  int nMc;
  float eMc[200],pMc[200];
  int idMc[200],statusMc[200];
    c.SetBranchAddress("nMc",&nMc);
    c.SetBranchAddress("energyMc",eMc);
    c.SetBranchAddress("pMc",pMc);
    c.SetBranchAddress("idMc",idMc);
    c.SetBranchAddress("statusMc",statusMc);

    c.SetBranchStatus("*",0);
    c.SetBranchStatus("*Mc",1);


  TFile outputFile(output.c_str(),"RECREATE");
  
  Long64_t iEntry=-1;
  Long64_t nEntries = -1;
  while(c.GetEntry(++iEntry)) {
    if(iEntry % 1000 ==0) std::cout << "Processing " << iEntry << " / " << nEntries << "\r" << std::flush;
    float m22=0,m23=0;
    for(int i=0;i<nMc;i++) {
      if(statusMc[i]==3 && idMc[i]==1000022) m22 = TMath::Sqrt(eMc[i]*eMc[i]-pMc[i]*pMc[i]);
      if(statusMc[i]==3 && idMc[i]==1000023) m23 = TMath::Sqrt(eMc[i]*eMc[i]-pMc[i]*pMc[i]);
      if(m22>0.01 && m23 >0.01) break;
    }
    N.Fill(m23,m22);
  }
  
  outputFile.cd();
  N.Write();
  outputFile.Close();

}
