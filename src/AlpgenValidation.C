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

// #define Vecbos
#include "ThiagoAnalysis.hh"
#include "Vecbos.hh"

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[150];
  char settingFileName[150];
  char outFileName[150];

  if ( argc < 3 ){
    cout << "Error at Input: please specify an input file including the list of input ROOT files" << endl; 
    cout << "                and an input file including the initialization values of the analysis parameters" << endl; 
    cout << "Example:        ./AlpgenValidation list.txt settings.txt output.root" << endl;
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  strcpy(settingFileName,argv[2]);
  strcpy(outFileName,argv[3]);

  TChain *theChain = new TChain("ntp1");
  char Buffer[500];
  char MyRootFile[2000];  
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	theChain->Add(MyRootFile);
	std::cout << "chaining " << MyRootFile << std::endl;
      }
  }
  inputFile->close();
  delete inputFile;

  theChain->MakeClass("thisiswhyitcrashed");

  ThiagoAnalysis vecbos(theChain);  
  vecbos.ReadParameters(settingFileName);
  vecbos.Loop(outFileName);
  system("\\rm thisiswhyitcrashed.h");
  system("\\rm thisiswhyitcrashed.C");

  return 0;
}
