//read the Hgg Energy Scale as a csv

// this class can do both energy scales and energy smears

#ifndef HggEnergyScale_hh
#define HggEnergyScale_hh

#include "ReadConfig.hh"
#include <string>
#include <map>
#include <vector>
#include "VecbosEGObject.hh"
//#include "TString.h"


class HggEnergyScale{
public:
  HggEnergyScale(std::string path);
  std::pair<float,float> getDEoE(VecbosPho,int);
  float getMCScaleErr(VecbosPho,int);
  bool isValid(){return valid;}
  int nRegions;
  std::vector<std::string> configNames;

  static const float r9Cut = 0.94;
  std::vector<bool> highR9;
  std::vector<float> minEta;
  std::vector<float> maxEta;
private:
  bool valid;
  
  float getCategory(VecbosPho);

  //energy scale
  std::vector<int> runs;
  std::map<std::string, std::vector<float> > energyScales; 
  
  int getRunIndex(int);
};
#endif

