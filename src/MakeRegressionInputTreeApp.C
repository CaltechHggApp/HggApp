//-------------------------------------------------------
// Description:
//    Routine to run Vecbos selection
// Authors:
//   Alex Mott (Caltech)
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>

// Vecbos includes
#include <CommonTools/include/TriggerMask.hh>
#include <include/Vecbos.hh>

#include "ArgParser.hh"
#include <include/MakeRegressionInputTree.hh>
#include <map>
#include <string>
using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  ArgParser a(argc,argv);

  a.addArgument("InputFile",ArgParser::required, "the input root file or list of input root files to process");
  a.addArgument("OutputFile",ArgParser::required, "the name of the file to write the output to");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  vector<string> fileNames;
  char Buffer[500];
  TString RootFileName;

  string inputFileName = a.getArgument("InputFile");

  if(inputFileName.find(".root") == string::npos) {

    ifstream *inputFile = new ifstream(a.getArgument("InputFile").c_str());
    vector<string> filesToRemove;
    while( !(inputFile->eof()) ){
      inputFile->getline(Buffer,2000);
      RootFileName = TString(Buffer).Strip();
      
      if(RootFileName.First('#') != -1 || RootFileName.Length()==0) continue;
      
    fileNames.push_back(string(RootFileName.Data()));
    
    std::cout << "chaining " << RootFileName << std::endl;
    } 
    inputFile->close();
    delete inputFile;
    
  }else{ // is a root file
    fileNames.push_back(inputFileName);
  }
  //  theChain->MakeClass("thisiswhyitcrashed");



  MakeRegressionInputTree sel(fileNames,"HggReduce", a.getArgument("OutputFile"));
  sel.setConfig("hgg2013/dummy.cfg");

  //sel.addTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v4");
  sel.Loop();
  
  cout << "DONE" <<endl;
  return 0;

}
