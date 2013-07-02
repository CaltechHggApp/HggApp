#ifndef MakeSpinToyWorkspace_h
#define MakeSpinToyWorkspace_h

#include "TObjArray.h"
#include "RooProdPdf.h"

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooGlobalFunc.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include "RooStats/SPlot.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooBernstein.h"
#include "RooLinkedListIter.h"
#include "RooCBShape.h"

#include "MakeSpinSeparationTest.h"
#include "MakeSpinPlots.h"
#include "MakeSpinFits.h"

class MakeSpinToyWorkspace : public MakeSpinSeparationTest{
public:
  MakeSpinToyWorkspace(TString fileName,TString wsName="cms_hgg_spin_workspace");

  bool saveWorkspaces;

  int Nbkg,Nsig;
  float lumi;

  TString tag;

  RooRealVar* mass;

  RooAbsPdf *bkgMassPdf,*hggMassPdf,*rsgMassPdf;

  void setup(TString t,float l);

  void generateN(int Ntoy);
  void generateToyWorkspace(bool doRSG);

  TObjArray toyWorkspaces;

  float getExpEventsCosT(float lumi, TString MC, TString cat,float minCosT, float maxCosT);

  double prob(float *exp, float *obs, float *err, int N);

  void save(TString filename);
};

#endif
