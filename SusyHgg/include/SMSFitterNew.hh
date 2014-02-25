#ifndef SMSFitterNew_hh
#define SMSFitterNew_hh

#include "FitterNew.hpp"

class SMSFitter : public Fitter {
public:
  SMSFitter(TString inputFileName,TString outputFileName);


  virtual ~SMSFitter();

  void setNEntriesFile(TString fileName);
  
  virtual void Run();

  static void getSMSPoints(std::vector<TString>* pts);

private:
  virtual void buildHistograms() override;
  virtual void processEntry() override;

  TString getSMSPoint();

  TFile* nEntriesFile=0;
  TH2F* nEntriesHist=0;

  std::vector<TString> sms_points;
};

#endif
