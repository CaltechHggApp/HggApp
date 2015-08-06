#ifndef FitterNew_hpp
#define FitterNew_hpp

#include "TString.h"
#include "TH2F.h"
#include "TH3F.h"
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
  Fitter(TString inputFileName,TString outputFileName,bool useHT=false);
  virtual ~Fitter();

  std::map<TString,float> xsecs = {
    {"gg_H_125",19.27},{"wz_H_125",0.7046},{"vbf_H_125",1.578},{"sms",0.3}
  };

  enum kSelectionSet {kLoose, kAN239, kHighPt};

  void setSelection(kSelectionSet set) { basicSelection = set; }

  void Run();

  bool passBasicSelection(); //!< event passes the basic event selection

  TString getCategory(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag,float mbbH, float mbbZ,float r9_1,float r9_2); //!< get the category for this event

  static TString getCategoryOrig(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag,float mbbH, float mbbZ,float r9_1,float r9_2); //!< get the category for this event
  static TString getCategoryAlt(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag,float mbbH, float mbbZ,float r9_1,float r9_2); //!< get the category for this event


  
  void setSigEff(TString cat, float se) { sigmaEffectives[cat] = se; }

  void setLumi(float l) {lumi=l;}

  virtual void processEntry(bool doSyst=true);

  void setXSec(float x){target_xsec=x;}
  void setNTotal(int n){N_total=n;}

  const static std::vector<TString>* getCatNames() {return &catNames;}
  const static std::vector<TString>* getSysNames() {return &systematicNames;}

  void setNSigEffs(float n){nSigEffSignalRegion=n;}

  //static constexpr float minMgg = 100;
  static constexpr float minMgg = 103;
  static constexpr float maxMgg = 160;

  void setUseHT(bool b=true){useHT=b;}

  void setMetPhiSF(TString s) { MetPhiSF_file=s; }

  void setDoAlternateAnalysis(bool b=true){ alternateAnalysis=b; }

protected:
  TFile *outputFile;

  bool isSMS=true;

  kSelectionSet basicSelection = kLoose;

  //do the blessed alternate analysis
  bool alternateAnalysis= false;

  float hggSigStrength = 1.00;

  float lumi = 1;
  //bool processData = true;

  float weight = 1.; //weight for this event

  const float triggerEff=0.81;

  int N_total=1;  //total number of events in the sample (for normalization)
  float target_xsec=1; //target x-sec in pb
  const float HggBR = 2.28E-3; //Br H-->photons

  float nSigEffSignalRegion=2;

  TString MetPhiSF_file="";

  const static std::vector<TString> catNames;
  const static std::vector<TString> systematicNames;
  std::vector<TString> systematicDir   = {"Up","Down"};

  //histograms
  virtual void buildHistograms();
  std::map<TString, TH2F*> SignalRegionHistograms;
  std::map<TString, TH2F*> SignalRegionHistogramsFineBin;
  std::map<TString, TH2F*> SignalRegionHistogramsHTMET;
  std::map<TString, TH2F*> SignalRegionHistogramsFineBinHTMET;
  std::map<TString, TH1D*> mgg_dists;

  std::map<TString, TH3F*> SignalRegions3D;
  
  //sigma effectives per category
  std::map<TString, float> sigmaEffectives;

  //ranges for the binning

  //   const static int nXbins =5;
  //   double xBins[nXbins] = {0,200,400,1000,2000};
  //   const static int nYbins =6;
  //   double yBins[nYbins] = {0,0.05,0.1,0.2,0.5,1.0};

  bool useHT=false;

  const static int nXbins =23;
  double xBins[nXbins] = {150,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000};
  const static int nYbins =21;
  double yBins[nYbins] = {0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
  const static int nZbins =200;
  double zBins[nZbins];

  //define systematics
  virtual float getSysErrPho(float eta,float r9);
  std::map<TString,float> smearSys = {
    {"EBLow_lowR9",0.0002},{"EBLow_highR9",0.0002},{"EBHigh_lowR9",0.0003},{"EBHigh_highR9",0.0014},
    {"EELow_lowR9",0.0003},{"EELow_highR9",0.0007},{"EEHigh_lowR9",0.0004},{"EEHigh_highR9",0.0003}
  };

  std::vector<float> SFbErr_ptMax = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};

  std::vector<float> SFbErr_CSVL = {
      0.033299,
      0.0146768,
      0.013803,
      0.0170145,
      0.0166976,
      0.0137879,
      0.0149072,
      0.0153068,
      0.0133077,
      0.0123737,
      0.0157152,
      0.0175161,
      0.0209241,
      0.0278605,
      0.0346928,
      0.0350099 };

  std::vector<float> SFbErr_CSVM = {
    0.0415707,
    0.0204209,
    0.0223227,
    0.0206655,
    0.0199325,
    0.0174121,
    0.0202332,
    0.0182446,
    0.0159777,
    0.0218531,
    0.0204688,
    0.0265191,
    0.0313175,
    0.0415417,
    0.0740446,
    0.0596716 
  };
  std::map<TString,int> nSignal;
  std::map<TString,int> nTotal;

  float getBTagSF(float errHigh=0, float errSec=0);

  float getSFb(float pt, bool CSVL);
};

#endif
