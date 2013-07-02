//-------------------------------------------------------
// Description:
//    Routine to get total number of events in  achain
// Authors:
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

using namespace std;

int main(int argc, char* argv[]) {

  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[150];

  if ( argc < 2 ){
    cout << "Error at Input: please specify an input file including the list of input ROOT files" << endl; 
    return 1;
  }

  strcpy(inputFileName,argv[1]);

  TChain *theChain = new TChain("ntp1");
  char Buffer[500];
  char MyRootFile[2000];  
  cout << "input: " << inputFileName << endl;
  ifstream *inputFile = new ifstream(inputFileName);
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	theChain->Add(MyRootFile);
	cout << "chaining " << MyRootFile << endl;
      }
  }
  inputFile->close();
  delete inputFile;

  cout << "Number of events found: " << theChain->GetEntries() << endl;

  return 0;
}
