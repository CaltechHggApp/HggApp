#ifndef HggMCWeight_hh
#define HggMCWeight_hh

#include "string.h"
#include "ReadConfig.hh"
#include "TH1F.h"
#include "TFile.h"
#include <cmath>

class HggMCWeight{
public:
  HggMCWeight();
  void init(ReadConfig&);
  void setReweighting(std::string,bool); // this will fail an assertion if the reweighting isn't defined
  float getWeight();
  void setHiggsMass(float m){higgsMass=round(m);}
  void setPileup(float p){thisPileup=p;}
  static float round(float f){ return floor(f+0.5); }

private:
  //BranchingRatioFile
  ReadConfig branching;
  //XSec
  ReadConfig xsec;

  // define which reweighting to do
  typedef std::map<std::string,bool> flagmap_t;
  flagmap_t flags;
  // map of the function to call
  typedef float (HggMCWeight::*weightMethod_t)();
  typedef std::map<std::string, weightMethod_t> funcmap_t;
  funcmap_t functions;


  //define the functions to get the actual scales;
  float getXSec();
  float getBranching();
  float getPileup();

  float higgsMass;
  float thisPileup;
  TFile * puWeightFile;
  TH1F* pileupWeight;
};

#endif
