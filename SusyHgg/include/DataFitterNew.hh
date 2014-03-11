#ifndef DataFitterNew_hh
#define DataFitterNew_hh

#include "FitterNew.hpp"

#include "RooWorkspace.h"

class DataFitter : public Fitter {
public:
  DataFitter(TString inputFileName, TString outputFileName): Fitter(inputFileName,outputFileName) {
    isSMS=false;
  }

  static RealVar doFitGetScale(TTree *data,float sigEff,RooWorkspace* ws=0); //!< for real data, do the background fit and use that to compute the sideband->signal scale factor

  virtual void Run();

protected:

  void buildSidebandHistograms();

  virtual void processEntrySidebands();

  virtual float getSysErrPho(float eta,float r9);

  std::map<TString,RealVar> scales;

  std::map<TString,TH2F*> SidebandRegionHistograms;
  std::map<TString,TH2F*> SidebandRegionHistogramsFineBin;

  std::vector<RooWorkspace*> mggFitWorkspaces;

  std::map<TString,float> scaleSys = {
    {"EBLow_lowR9",0.0004},{"EBLow_highR9",0.0003},{"EBHigh_lowR9",0.0008},{"EBHigh_highR9",0.0008},
    {"EELow_lowR9",0.0013},{"EELow_highR9",0.0011},{"EEHigh_lowR9",0.0010},{"EEHigh_highR9",0.0009},
  };

};

#endif
