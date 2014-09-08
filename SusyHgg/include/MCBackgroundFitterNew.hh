#ifndef MCBackgroundFitterNew_hh
#define MCBackgroundFitterNew_hh

#include "DataFitterNew.hh"

#include "RooWorkspace.h"

class MCBackgroundFitter : public DataFitter {
public:
  MCBackgroundFitter(TString inputFileName, TString outputFileName,bool useHT): DataFitter(inputFileName,outputFileName,useHT) {
    isSMS=false;
  }

  virtual void Run() override;

  virtual void fixNorm(TString catName, float norm) override;

  void setCorrSherpa(bool b=true){correctSherpaEnhance=b;}
protected:

  virtual void buildSidebandHistograms() override;

  void calcNormWeights();
  std::map<TString,float> normMap; //the number of events we want in the signal region of each cat
  std::map<TString,float> weightMap; //the scale factor needed to get that # events

  bool correctSherpaEnhance=false;

  float enhance_3jet=2.0;
  float enhance_4jet=5.0;
  float getSherpaCorrection();
};

#endif

