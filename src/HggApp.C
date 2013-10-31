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

#include <include/HggReducer.hh>
#include <include/HggMakePhotonTree.hh>

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char** argv) {
  /// Gets the list of input files and chains
  /// them into a single TChain
  ArgParser a(argc,argv);

  a.addArgument("InputList",ArgParser::required,"list of input (vecbos) nTuples");
  a.addArgument("OutputFile",ArgParser::required,"path to the output file");
  a.addArgument("ConfigFile",ArgParser::required,"input configuration file");

  a.addLongOption("isData",ArgParser::noArg,"specify that this is data (rather than MC)");
  a.addLongOption("json",ArgParser::reqArg,"specify the JSON for data");
  a.addLongOption("start",ArgParser::reqArg,"specify the start position in the tree");
  a.addLongOption("stop",ArgParser::reqArg,"specify the stop position in the tree");

  a.addLongOption("XRootD",ArgParser::reqArg,"access input files with XRootD instead of the specified path [argument is the site]");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }
  
  TChain *theChain = new TChain("ntp1");
  char Buffer[2000];
  TString RootFileName;
  char MyRootFile[2000];  
  ifstream *inputFile = new ifstream(a.getArgument("InputList").c_str());
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
  //delete inputFile;
  // get additional input options
  int start = 0; 
  int stop  = -1;
  bool isData = a.longFlagPres("isData");
  string json = "";
  if(a.longFlagPres("start")) start = atoi(a.getLongFlag("start").c_str());
  if(a.longFlagPres("stop")) stop = atoi(a.getLongFlag("stop").c_str());
  if(a.longFlagPres("json")) json = a.getLongFlag("json");
  
  cout << "Running on: " << stop-start << " Entries" << endl;

  HggReducer vecbos(theChain, json, isData, isData);
  cout << "Setting Config File:" << endl;
  cout << string(a.getArgument("ConfigFile")) << endl;
  vecbos.setConfig(a.getArgument("ConfigFile"));
  vecbos.Loop(a.getArgument("OutputFile"), start, stop);

  return 0;
}
