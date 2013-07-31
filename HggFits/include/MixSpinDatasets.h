#ifndef MixSpinDatasets_h
#define MixSpinDatasets_h

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooCategory.h"

#include "TString.h"

#include <vector>
#include <string>
#include "MakeSpinFits.h"

class MixSpinDatasets{
public:
  MixSpinDatasets(RooWorkspace *w);

  void scheduleMix(const char* mc1, const char* mc2, float f1,TString outputName="");

  void scheduleMerge(std::vector<std::string> mcNames, TString outputName="");

  void mixAll();
  void mergeAll();

  void mix(const char* mc1, const char* mc2, float f1,TString outputName="");
  void merge(std::vector<TString> names,TString outputName="");
protected:
  RooWorkspace *ws;

  std::vector<TString> catNames;

  std::vector<TString> mc1L;
  std::vector<TString> mc2L;
  std::vector<float> f1L;
  std::vector<TString> outputNameL;

  std::vector<std::vector<TString> > mergeL;
  std::vector<TString> mergeOutputNameL;

  void internalMix(const char* mc1, const char* mc2, float f1,TString outputName,TString cat);
  void internalMerge(std::vector<TString> names, TString outputName,TString cat,RooCategory* labels);
};

#endif
