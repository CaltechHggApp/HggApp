#ifndef Fitter_hpp
#define Fitter_hpp

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

#include "TString.h"
#include "TH2F.h"
#include "TFile.h"

#include <vector>
#include <map>
#include <exception>
#include <stdexcept>

struct fitInfo {
  float sideband_low_min;
  float sideband_low_max;

  float signal_min;
  float signal_max;

  float sideband_high_min;
  float sideband_high_max;

};

class Fitter {
public:
  Fitter(TString inputFileName,TString outputFileName);
  ~Fitter();

  std::map<TString,float> xsecs = {
    {"gg_H_125",19.27},{"wz_H_125",0.7046},{"vbf_H_125",1.578},{"sms",0.3}
  };

  virtual void Run();

  virtual void buildHistograms(int catIndex);
  virtual void runFits(int catIndex);
  virtual void computeScaleFactor(int catIndex);
  virtual void buildSubtractedHistograms(int catIndex);

  TH2F* getSignalRegionHistogram(const TH2F& hist,const char* name);

  void addMcName(TString name);

private:
  TFile *inputFile;
  TFile *outputFile;
  RooWorkspace *ws;

  std::vector<TString> cats;
  std::vector<fitInfo> per_cat_fit_ranges;
  std::vector<std::pair<float,float>> sideband_integrals;
  std::vector<std::pair<float,float>> signal_integrals;
  std::vector<std::pair<float,float>> sideband_to_signal_scale_factors;

  std::vector<TString> sm_higgs_names;
  std::vector<TString> sms_names;
  TString data_name="data";

  enum region{kSignal,kSideband,kAll};
  TH2F* build_histogram(TString name, int catIndex,region reg,float weight=1);

  std::map<TString,TH2F*> save_histograms;
  void saveAll();


  void defineCats();

  virtual void setupCategory(TString name);
};

#endif
