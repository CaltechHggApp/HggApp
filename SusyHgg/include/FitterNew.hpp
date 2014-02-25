#ifndef FitterNew_hpp
#define FitterNew_hpp

#include "TString.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "SusyHggTree.h"

#include <vector>
#include <map>
#include <exception>
#include <stdexcept>

#include "assert.h"

struct fitInfo {
  float sideband_low_min;
  float sideband_low_max;

  float signal_min;
  float signal_max;

  float sideband_high_min;
  float sideband_high_max;

};

struct RealVar {
  float val;
  float error;
};

/*
  takes in a TTree 
*/

class Fitter : public SusyHggTree{
public:
  Fitter(TString inputFileName,TString outputFileName);
  virtual ~Fitter();

  std::map<TString,float> xsecs = {
    {"gg_H_125",19.27},{"wz_H_125",0.7046},{"vbf_H_125",1.578},{"sms",0.3}
  };

  void Run();

  bool passBasicSelection(); //!< event passes the basic event selection


  TString getCategory(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag); //!< get the category for this event

  
  void setSigEff(TString cat, float se) { sigmaEffectives[cat] = se; }

  void setLumi(float l) {lumi=l;}

  virtual void processEntry();

  void setXSec(float x){target_xsec=x;}
  void setNTotal(int n){N_total=n;}

  const static std::vector<TString>* getCatNames() {return &catNames;}
  const static std::vector<TString>* getSysNames() {return &systematicNames;}


protected:
  TFile *outputFile;

  bool isSMS=true;

  float lumi = 1;
  //bool processData = true;

  float weight = 1.; //weight for this event

  int N_total=1;  //total number of events in the sample (for normalization)
  float target_xsec=1; //target x-sec in pb
  const float HggBR = 2.28E-3; //Br H-->photons


  const static std::vector<TString> catNames;
  const static std::vector<TString> systematicNames;
  std::vector<TString> systematicDir   = {"Up","Down"};

  //histograms
  virtual void buildHistograms();
  std::map<TString, TH2F*> SignalRegionHistograms;
  std::map<TString, TH2F*> SignalRegionHistogramsFineBin;
  std::map<TString, TH1D*> mgg_dists;
  
  //sigma effectives per category
  std::map<TString, float> sigmaEffectives;

  //ranges for the binning
  const static int nXbins =5;
  double xBins[nXbins] = {0,200,400,1000,2000};
  const static int nYbins =6;
  double yBins[nYbins] = {0,0.05,0.1,0.2,0.5,1.0};

  //define systematics
  virtual float getSysErrPho(float eta,float r9);
  std::map<TString,float> smearSys = {
    {"EBLow_lowR9",0.0002},{"EBLow_highR9",0.0002},{"EBHigh_lowR9",0.0003},{"EBHigh_highR9",0.0014},
    {"EELow_lowR9",0.0003},{"EELow_highR9",0.0007},{"EEHigh_lowR9",0.0004},{"EEHigh_highR9",0.0003}
  };

  std::map<float,float> SFb_error = { //min_Pt --> error
    {20,0.0484285},{30,0.0126178},{40,0.0120027},{50,0.0120027},{60,0.0145441},{60,0.0131145},{70,0.0168479},
    {80,0.0160836},{100,0.0126209},{120,0.0136017},{160,0.019182},{210,0.0198805},{320,0.0386531},{400,0.0392831},
    {500,0.0481008},{600,0.0474291}
  };


};

#endif
