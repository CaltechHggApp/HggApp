//-------------------------------------------------------
// Description:
//    Routine to run Vecbos selection
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//    Maurizio Pierini
//    CERN
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

#include <include/HggEfficiencyMap.hh>
#include <map>
using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[400];
  char outFileName[400];
  char cfg[400];
  if ( argc < 3 ){
    cout << "Error at Input: please specify an input file including the list of input ROOT files" << endl; 
    cout << "Example:        ./VecbosApp list.txt output.root cfg" << endl;
    cout << "Example:        ./VecbosApp list.root output.root cfg" << endl;

    return 1;
  }

  std::map<std::string,bool> Commands;
  Commands["--mmg"]=false;
  Commands["--zeeMC"]=false;
  Commands["--zee"]=false;

  for(int arg=3; arg<argc; arg++){
    for(std::map<std::string,bool>::iterator it=Commands.begin(); it!=Commands.end(); it++){
      if(it->first.compare(argv[arg])==0) it->second=true;
    }
  }

  for(std::map<std::string,bool>::iterator it=Commands.begin(); it!=Commands.end(); it++){
    if(it->second) cout << "Doing " << it->first << endl;
  }  
  // rad running options
  strcpy(inputFileName,argv[1]);
  strcpy(outFileName,argv[2]);
  strcpy(cfg,argv[3]);


  vector<string> fileNames;
  char Buffer[500];
  char MyRootFile[2000];  
  ifstream *inputFile = new ifstream(inputFileName);
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
  // get additional input options
  HggEfficiencyMap sel(fileNames,"HggReduce", string(outFileName));
  //sel.addTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v4");
  sel.setConfigFile(cfg);
  if(Commands["--mmg"]) sel.setMuMuGamma();
  if(Commands["--zeeMC"]) sel.setZeeMC();
  if(Commands["--zee"]) sel.setZeeData();
  sel.Loop();
  
  cout << "DONE" <<endl;
  system("rm thisiswhyitcrashed*");

  return 0;

}
