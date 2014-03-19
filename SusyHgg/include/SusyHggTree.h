#ifndef SusyHggTree_h
#define SusyHggTree_h

#include "SusyHggTreeBase.h"


class SusyHggTree : public SusyHggTreeBase {
public :
  SusyHggTree(TString fileName);
};

SusyHggTree::SusyHggTree(TString fileName):SusyHggTreeBase() {
  TFile *inputFile = new TFile(fileName);
  TTree *tree = (TTree*)inputFile->Get("SusyHggTree");
  Init(tree);
}

#endif
