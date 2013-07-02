#ifndef HggScaling_hh
#define HggScaling_hh

#include <VecbosEGObject.hh>
#include <ReadConfig.hh>
#include <string>
#include <vector>
#include <map>

class HggScaling{
public:
  HggScaling();
  HggScaling(std::string);
  void LoadConfig(std::string);
  void ScalePhoton(VecbosPho&);
private:
  // map between pair(Correction Name, IsEndcap) and pair(Linear,Constant) correction
  typedef std::map<std::pair<std::string,bool>,std::pair<float,float> > CorrectionMap;
  CorrectionMap Corrections;
  

  void ApplyScaling(VecbosPho&,std::pair<std::string,bool>);
};


#endif
