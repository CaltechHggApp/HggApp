#ifndef VecbosJetID_hh
#define VecbosJetID_hh

#include "VecbosEGObject.hh"
#include <map>

class VecbosJetID {
public:
  VecbosJetID();

  enum WP{kLoose,kMedium,kTight};

  bool passID(const VecbosJet& jet,WP wp);
  
protected:
  std::map<VecbosJetID::WP,float> eFracCut;

};

#endif
