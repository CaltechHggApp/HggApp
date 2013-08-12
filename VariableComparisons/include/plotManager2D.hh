#ifndef plotManager2D_hh
#define plotManager2D_hh

#include "TH2F.h"
#include "plotManager.hh"

class plotManager2D: public plotManager{
public:
  plotManager2D(TString inTag=""): plotManager(inTag){}
  std::array<TH2F*,3> get2DHistogram(TString cat, TString var1,TString var2);
  virtual void processChain(TChain *fChain,float weight);
  void saveAll2D(TFile *f);

protected:
  std::vector<std::vector<std::vector<TH2F*>>> fake2DHistograms;
  std::vector<std::vector<std::vector<TH2F*>>> realPho2DHistograms;
  std::vector<std::vector<std::vector<TH2F*>>> realEle2DHistograms;

  void buildHistograms2D();
  void processEntry2D(float weight);

};

#endif
