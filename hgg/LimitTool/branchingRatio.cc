#ifndef branchingRatio_cc
#define branchingRatio_cc
  //copied from the twiki: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
#include <TString.h>
class branchingRatio{
public:
  branchingRatio();
  std::map<float,float> BranchingRatios;
  std::map<TString,std::map<float,float> >Proc8TeV;
  std::map<float,float> Proc8TeV_GluGlu;
  std::map<float,float> Proc8TeV_VBF;
};
#endif
branchingRatio::branchingRatio(){
  BranchingRatios[90.]  = 1.23E-03;
  BranchingRatios[95.]  = 1.40E-03;
  BranchingRatios[100.] = 1.59E-03;
  BranchingRatios[105.] = 1.78E-03;
  BranchingRatios[110.] = 1.97E-03;
  BranchingRatios[115.] = 2.13E-03;
  BranchingRatios[120.] = 2.25E-03;
  BranchingRatios[123.] = 2.28E-03;
  BranchingRatios[124.] = 2.29E-03;
  BranchingRatios[125.] = 2.29E-03;
  BranchingRatios[129.] = 2.27E-03;
  BranchingRatios[130.] = 2.26E-03;
  BranchingRatios[135.] = 2.13E-03;
  BranchingRatios[140.] = 1.93E-03;
  BranchingRatios[145.] = 1.67E-03;
  BranchingRatios[150.] = 1.36E-03;
  BranchingRatios[155.] = 0.999E-03;

  Proc8TeV_GluGlu[90.]  = 36.772;
  Proc8TeV_GluGlu[95.]  = 33.153;
  Proc8TeV_GluGlu[100.] = 30.089;
  Proc8TeV_GluGlu[105.] = 27.360;
  Proc8TeV_GluGlu[110.] = 25.012;
  Proc8TeV_GluGlu[115.] = 22.937;
  Proc8TeV_GluGlu[120.] = 21.109;
  Proc8TeV_GluGlu[123.] = 20.113;
  Proc8TeV_GluGlu[124.] = 19.796;
  Proc8TeV_GluGlu[125.] = 19.487;
  Proc8TeV_GluGlu[129.] = 18.317;
  Proc8TeV_GluGlu[130.] = 18.040;
  Proc8TeV_GluGlu[135.] = 16.745;
  Proc8TeV_GluGlu[140.] = 15.578;
  Proc8TeV_GluGlu[145.] = 14.525;
  Proc8TeV_GluGlu[150.] = 13.567;
  Proc8TeV_GluGlu[155.] = 12.678;

  Proc8TeV_VBF[100] = 1.971;
  Proc8TeV_VBF[110] = 1.791;
  Proc8TeV_VBF[115] = 1.709;
  Proc8TeV_VBF[120] = 1.632;
  Proc8TeV_VBF[123] = 1.588;
  Proc8TeV_VBF[125] = 1.559;
  Proc8TeV_VBF[130] = 1.490;
  Proc8TeV_VBF[135] = 1.425;
  Proc8TeV_VBF[140] = 1.365;
  Proc8TeV_VBF[145] = 1.306;
  Proc8TeV_VBF[150] = 1.251;
 
  Proc8TeV["ggh"] = Proc8TeV_GluGlu;
  Proc8TeV["vbf"] = Proc8TeV_VBF;

}
