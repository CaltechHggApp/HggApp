#ifndef MixSpinDatasets_h
#define MixSpinDatasets_h

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooCategory.h"

#include "TString.h"

#include <vector>
#include "MakeSpinFits.h"

class MixSpinDatasets{
public:
  MixSpinDatasets(RooWorkspace *w);

  void scheduleMix(const char* mc1, const char* mc2, float f1,TString outputName="");

  void mixAll();

  void mix(const char* mc1, const char* mc2, float f1,TString outputName="");
protected:
  RooWorkspace *ws;

  std::vector<TString> catNames;

  std::vector<TString> mc1L;
  std::vector<TString> mc2L;
  std::vector<float> f1L;
  std::vector<TString> outputNameL;

  void internalMix(const char* mc1, const char* mc2, float f1,TString outputName,TString cat);
};

#endif
