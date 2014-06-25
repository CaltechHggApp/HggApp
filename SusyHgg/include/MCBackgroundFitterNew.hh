#ifndef MCBackgroundFitterNew_hh
#define MCBackgroundFitterNew_hh

#include "DataFitterNew.hh"

#include "RooWorkspace.h"

class MCBackgroundFitter : public DataFitter {
public:
  MCBackgroundFitter(TString inputFileName, TString outputFileName): DataFitter(inputFileName,outputFileName) {
    isSMS=false;
  }

  virtual void Run() override;

  virtual void fixNorm(TString catName, float norm) override;

protected:

  virtual void buildSidebandHistograms() override;

  void calcNormWeights();
  std::map<TString,float> normMap; //the number of events we want in the signal region of each cat
  std::map<TString,float> weightMap; //the scale factor needed to get that # events
};

#endif

