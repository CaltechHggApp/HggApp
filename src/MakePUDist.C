//
//  Simple Utility to make the PU distribution from reduced or vecbos nTuples
//
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"

#include "ArgParser.hh"

int main(int argc, char** argv) {

  ArgParser a(argc,argv);
  
  a.addArgument("InputList",ArgParser::required,"list of input nTuples");
  a.addArgument("OutputFile",ArgParser::required,"output file");
  a.addArgument("TreeName",ArgParser::optional,"optional tree name [ntp1]");

  std::string ret;
  if(a.process(ret) !=0){
    std::cout << "Invalid Options:  " << ret << std::endl;
    a.printOptions(argv[0]);
    return 0;
  }


  std::string treeName = a.getArgument("");
  if(treeName.compare("")==0) treeName = "ntp1";

  TChain *theChain = new TChain(treeName.c_str());
  char Buffer[2000];
  char MyRootFile[2000];  
  std::ifstream *inputFile = new std::ifstream(a.getArgument("InputList").c_str());
  if(!inputFile){
    std::cout << "Invalid input file" <<std::endl;
    return 1;
  }
  // get the tree with the conditions from the first file
  //  TTree *treeCond = new TTree();
  //  int nfiles=1;
  char tmpFileName[2000];
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,2000);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	if(std::string(MyRootFile).find("eos") != std::string::npos) {
	  theChain->Add(TString(MyRootFile));
        } else if(std::string(MyRootFile).find("castor") != std::string::npos) {
	  theChain->Add("rfio:"+TString(MyRootFile));
	} else{
	  theChain->Add(TString(MyRootFile));	 
	}
	std::cout << "chaining " << MyRootFile << std::endl;
      }
  } 
  inputFile->close();

  TFile *outputFile = new TFile(a.getArgument("OutputFile").c_str(),"RECREATE");
  TH1D  *outputHist = new TH1D("pileup_hist","",1600,0,160);
  Long64_t iEntry=-1;

  theChain->SetBranchStatus("*",0);
  theChain->SetBranchStatus("nPU",1);
  int nPU[3];
  theChain->SetBranchAddress("nPU",nPU);
  
  while(theChain->GetEntry(++iEntry)) {
    outputHist->Fill(nPU[1]);
  }
  outputHist->Write();
  outputFile->Close();
    

}
