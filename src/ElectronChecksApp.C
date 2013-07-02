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

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "include/ElectronChecks.hh"

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files

int main(int argc, char* argv[]) 
{
  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[150];

  char outFileName[150];

  if ( argc < 3 )
    {
      cout << "Error at Input: please specify an input file including the list of input ROOT files" << endl; 
      cout << "                the prefix for the output " << endl; 
      cout << "Example:        ./ElectronChecks list.txt eleChecks" << endl;
      return 1;
    }
  
  strcpy(inputFileName,argv[1]);
  strcpy(outFileName,argv[2]);
  
  TChain *theChain = new TChain("ntp1");
  char Buffer[500];
  char MyRootFile[2000];  
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);

  // get the tree with the conditions from the first file
  TTree *treeCond = new TTree();
  int nfiles=1;

  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	theChain->Add(MyRootFile);
        if ( nfiles==1 ) {
          TFile *firstfile = TFile::Open(MyRootFile);
          treeCond = (TTree*)firstfile->Get("Conditions");
        }
	std::cout << "chaining " << MyRootFile << std::endl;
        nfiles++;
      }
  }
  inputFile->close();
  delete inputFile;

  // for electron analysis:
  ElectronChecks eleChecks(theChain);
//   TriggerMask mask(treeCond);
//   mask.requireTrigger("HLT_Ele10_SW_L1R");
//   std::vector<int> requiredTriggers = mask.getBits();
//   vecbos.requireTrigger(requiredTriggers);
  eleChecks.Loop();
  eleChecks.writeAllHistos(TString(outFileName)+".root");
  
  return 0;
}
