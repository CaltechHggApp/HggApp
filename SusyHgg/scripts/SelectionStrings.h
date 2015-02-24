#ifndef SelectionStrings_h
#define SelectionStrings_h
#include "TString.h"

class SelectionStrings {
public:
  SelectionStrings();
  TString baseSelection = "(pho1_pt>40 && pho2_pt>25 && abs(pho1_eta)<1.48 && abs(pho2_eta)<1.48 && pho1_pass_iso && pho2_pass_iso)";
  TString mggSigRegion[5];

  TString boxDefs[5];
};

SelectionStrings::SelectionStrings() {
  mggSigRegion[0] = "(mgg>121.96 && mgg<129.04)";
  mggSigRegion[1] = "(mgg>121 && mgg<130)";
  mggSigRegion[2] = "(mgg>121 && mgg<130)";
  mggSigRegion[3] = "(mgg>122.04 && mgg<128.96)";
  mggSigRegion[4] = "(mgg>120 && mgg<131)";

  boxDefs[0] = "(ptgg>110)";
  boxDefs[1] = "(ptgg<=110 && abs(mbb_NearH-125)<15)";
  boxDefs[2] = "(ptgg<=110 && abs(mbb_NearH-125)>15 && abs(mbb_NearZ-91.2)<15)";
  boxDefs[3] = "(ptgg<=110 && abs(mbb_NearH-125)>15 && abs(mbb_NearZ-91.2)>15 && pho1_sigEoE<0.015 && pho2_sigEoE<0.015)";
  boxDefs[4] = "(ptgg<=110 && abs(mbb_NearH-125)>15 && abs(mbb_NearZ-91.2)>15 && !(pho1_sigEoE<0.015 && pho2_sigEoE<0.015))";

}

#endif
