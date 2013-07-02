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
#include <include/HggSelector.hh>
#include <map>
#include <string>
using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  ArgParser a(argc,argv);

  a.addArgument("InputFile",ArgParser::required, "the input root file or list of input root files to process");
  a.addArgument("OutputFile",ArgParser::required, "the name of the file to write the output to");
  a.addArgument("ConfigFile",ArgParser::required, "the name of the configuration file");

  a.addLongOption("isData",ArgParser::noArg,"specify that this is real data");
  a.addLongOption("doMuMuGamma",ArgParser::noArg,"Write mu mu gamma output");
  a.addLongOption("suppressElectronVeto",ArgParser::noArg,"don't veto photons matched to electrons (for Zee)");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  vector<string> fileNames;
  char Buffer[500];
  char MyRootFile[2000];  
  ifstream *inputFile = new ifstream(a.getArgument("InputFile").c_str());
  // get the tree with the conditions from the first file
  //  TTree *treeCond = new TTree();
  //  int nfiles=1;
  char tmpFileName[256];
  vector<string> filesToRemove;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	if(string(MyRootFile).find("eos") != std::string::npos) {
	  fileNames.push_back(string(MyRootFile));
        } else if(string(MyRootFile).find("castor") != std::string::npos) {
	  fileNames.push_back("rfio:"+string(MyRootFile));
	} else{
	  fileNames.push_back(string(MyRootFile));	 
	}
        // theChain->Add("root://castorcms/"+TString(MyRootFile));
	//        theChain->Add(TString(MyRootFile));
	std::cout << "chaining " << MyRootFile << std::endl;
	//	if ( nfiles==1 ) {
	//	  TFile *firstfile = TFile::Open("root://castorcms/"+TString(MyRootFile));
	//	  treeCond = (TTree*)firstfile->Get("Conditions");
	//	}
	//        nfiles++;
      }
  }
 

  //  theChain->MakeClass("thisiswhyitcrashed");

  inputFile->close();
  delete inputFile;


  HggSelector sel(fileNames,"HggReduce", a.getArgument("OutputFile"));
  //sel.addTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v4");
  sel.setConfig(a.getArgument("ConfigFile"));
  if(a.longFlagPres("doMuMuGamma")) sel.setDoMuMuGamma();
  if(a.longFlagPres("isData")) sel.setIsData(true);
  else sel.setIsData(false);
  if(a.longFlagPres("suppressElectronVeto")) sel.suppressElectronVeto();

  sel.Loop();
  
  cout << "DONE" <<endl;
  system("rm thisiswhyitcrashed*");

  return 0;

}
