/**

Creates a Signal Injection test for the new Fitter Code

*/

#ifndef SignalInjectionNew_hh
#define SignalInjectionNew_hh

#include "TH2F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TString.h"

#include <vector>

#include "FitterNew.hpp"
#include "TList.h"


class SignalInjection {
public:
  SignalInjection(TString outputFileName);
  virtual ~SignalInjection();

  void setDataFile(TString name) {dataFile=openInjFile(name,dataFile);}
  void setSMSFile(TString name) {smsFile=openInjFile(name,smsFile);}
  void addHiggsFile(TString name);

  void setInjXsec(float x){xsec=x;}
  void setMassPoint(TString m){massPoint=m;}


  void make();

  const TFile* getDataFile(){return dataFile;}

protected:
  TFile *outputFile=0;

  TFile * dataFile=0;
  TFile * smsFile=0;
  std::vector<TFile *> smHiggsFiles;
  
  TH2F* getToy(const TH2F& input,TString name="");

  TRandom3 rng;

  TFile* openInjFile(TString fileName, TFile *file);

  const std::vector<TString> *cats = Fitter::getCatNames();
  
  float xsec =1 ;
  TString massPoint="";
  

};

#endif
