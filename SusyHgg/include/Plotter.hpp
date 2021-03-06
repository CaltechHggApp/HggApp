#ifndef Plotter_hpp
#define Plotter_hpp

#include "RooWorkspace.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"
#include "RooBinning.h"
#include "RooGaussian.h"
#include "RooLognormal.h"

#include "TString.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TList.h"
#include "TLatex.h"

#include <vector>
#include <map>
#include <exception>
#include <stdexcept>
#include "assert.h"

class Plotter {
public:
  Plotter(TString inputFileName, TString outDir, TString outTag);
  ~Plotter();

  void plotTH2s();
  void plotMggFits(TString catName);
  void plotSignalPeaks(TString catName);
  void plotFrenchFlag(TString catName);

  void Run();

  void setIsSMS(bool b=true){isSMS=b;}

private:
  RooWorkspace* inputWs;
  TFile* inputFile;

  TString outputTag;
  TString outputDir;

  void getCategories(RooCategory* roocat);
  std::vector<TString> categories;

  RooRealVar *mgg;

  void setStyle();
  TStyle *vecbosStyle;

  bool isSMS=false;
};


#endif
