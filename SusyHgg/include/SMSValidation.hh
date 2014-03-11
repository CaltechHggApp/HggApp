/**

perform a set of validation exercises on the SMSs for FastSIM versus FullSIM

*/

#ifndef SMSValidation_hh
#define SMSValidation_hh

#include "TString.h"
#include "FitterNew.hpp"
#include "TH1.h"
#include "SusyHggTree.h"

#include <map>


class SMSValidation: public SusyHggTree {
public:
  SMSValidation(TString inputFileName, TString outputFileName);
  virtual ~SMSValidation();

  virtual void Run();

protected:
  TFile *outputFile;
  std::map<TString,TH1*> plots;

  virtual void setupPlots();
  virtual void processEntry();

  const std::vector<TString>* catNames = Fitter::getCatNames();
  
};

#endif
