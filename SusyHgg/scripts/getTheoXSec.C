#include "TGraphErrors.h"
#include "TString.h"
#include "TTree.h"

TGraphErrors* getTheoXSec(TString file) {
  TTree t;

  t.ReadFile(file);

  TGraphErrors *g = new TGraphErrors(t.GetEntries());

  float m,x,e;
  t.SetBranchAddress("mass",&m);
  t.SetBranchAddress("xsec",&x);
  t.SetBranchAddress("err",&e);
  

  int iEntry=-1;
  while(t.GetEntry(++iEntry)) {
    g->SetPoint(iEntry,m,x);
    g->SetPointError(iEntry,0,e);
  }
  return g;
}
