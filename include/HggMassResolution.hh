#ifndef HggMassResolution_hh
#define HggMassResolution_hh

#include <VecbosEGObject.hh>
#include <string>
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

struct diPhotonInfo{
  TLorentzVector p4Pho1;
  TLorentzVector p4Pho2;
  

  float eRes1;
  float eRes2;
  pair<float,float> catRes1;
  pair<float,float> catRes2; // resolutions for each photon category
  double diPhoMass;
  double angle;
};

class HggMassResolution{
public:
  HggMassResolution();
  void setConfigFile(string s){config = s;}
  void setCategoryType(string s){categoryType = s;}
  void init();
  double getMassResolution(VecbosPho*,VecbosPho*,TVector3,bool);
  double getMassResolutionEonly(VecbosPho*,VecbosPho*,TVector3);
  const static int nCategories=9;
  std::vector<string>Categories;
  const static float r9Cut=0.94;
  std::vector<bool> highR9;
  std::vector<float> minEta;
  std::vector<float> maxEta;
  std::vector<float> dzRes;
  const static int sphericalIndex=1;
  float getResolution(VecbosPho*);
private:
  void clear();
  string config;
  string categoryType;
  double diPhoMass;
  int getCategory(VecbosPho*);
  double getAngleResolution(VecbosPho*, VecbosPho*, TVector3, bool);
  bool isSphericalPhoton(int,int);

  pair<float,float> smear[nCategories];

};


#endif 
