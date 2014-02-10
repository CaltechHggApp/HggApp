#ifndef SimpleFitter_hpp
#define SimpleFitter_hpp

#include "Fitter.hpp"

class SimpleFitter : public Fitter {
public:
  SimpleFitter(TString inputFileName, TString outputFileName) : Fitter(inputFileName,outputFileName) {}

  virtual void Run();
  virtual void buildHistograms(int catIndex);
  virtual void SMFit();

protected:

  virtual TH2F* build_histogram(TString name, int catIndex, region reg, float weight=1);
}

#endif
