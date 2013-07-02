#ifndef CUTBASEDSELECTOREE_HH
#define CUTBASEDSELECTOREE_HH

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include <vector>

struct SelectorData {
  //! all the variables used to do the selection
  int signal;
  int inacceptance;
  int njetsMC, njets;
  bool passedHLT,matchedHLT;
  int nRecoEle, nAccEle, nIdTightEle, nIdLooseEle,
    nIsolTightEle, nIsolLooseEle, nConvRejTightEle, nConvRejLooseEle,
    nPV, nDzVertexEle, nDxyVertexEle, nMuons;
  float btagEVT, mInv, mhtMET;
  bool chargeMajorityMethod;
  
  //! W specific cuts
  float met, mt;
  bool foundAnyZ;
  std::vector<float> mhtJet;

  //! event weight
  float weight;

  //! contructor
  SelectorData();
  //! copy constructor
  SelectorData(SelectorData &data);
};

class CutBasedSelectorEE {

public:
  CutBasedSelectorEE(SelectorData data);
  ~CutBasedSelectorEE() {};
  
  //! set true if it is MC (by default is data)
  void isMc(int what=1) {mc_ = what;} ;
  //! output of the selection for Z
  bool outputZ(Selection* commonSelection, Selection* selection,
              Counters* inclusiveCounter,
              Counters (*mcJetBinCounter)[6], Counters (*recoJetBinCounter)[6]);
  //! output of the selection for W
  bool outputW(Selection* commonSelection, Selection* selection,
               Counters* inclusiveCounter,
               Counters (*mcJetBinCounter)[6], Counters (*recoJetBinCounter)[6]);

private:
  
  SelectorData data_;
  int mc_;
};

#endif
