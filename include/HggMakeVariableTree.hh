#include <VecbosEGObject.hh>

#include <vector>
#include <iostream>
#include <string>

#include <TChain.h>
#include "BaseSelector.hh"

#include "HggPhotonID.hh"

class HggMakeVariableTree : public BaseSelector{
public:
  HggMakeVariableTree(vector<string> fNames,string treeName,string outputFile):BaseSelector(fNames,treeName,outputFile){setDoFill(false);}

  inline void setMinPt(float p) {minPt = p;}
  inline void setRequireGenMatchPhoton(bool b) { requireGenMatchPhoton=b;}
  inline void setRequireGenMatchElectron(bool b) { requireGenMatchElectron=b;}
protected:
  //mandatory overrides
  virtual void processEntry(Long64_t iEntry);
  virtual int init(){ return 0; }
  virtual void processConfig(ReadConfig &cfg);
  virtual void setupOutputTree();
  
  virtual void clear();
  virtual void firstInit();

  HggPhotonID *photonID;

  //selection criteria
  float minPt = -1;

  //output variables
   float se;
   float etaSC;
   float r9;
   float phi;
   float pt;
   int realPho;
   int realEle;
   int passPre;
   int pos;

   float pfChargedGood;
   float pfChargedWorst;
   float pfPhoton;
   float HE;
   float sieie;
   float sieip;
   float sipip;

   int passCiC;
   int passCiC_id,passCiC_iso;
   int Trigger;
   int TightPt;
   float mass;

   float e3x3;
   float e5x5;
   float rawE;
   float etaWidth;
   float phiWidth;
   int nBC;
   float energyBC;
   float etaBC;
   float phiBC;

   float eMax;
   float e2nd;
   float eTop;
   float eBottom;
   float eLeft;
   float eRight;

   float e2x5Max;
   float e2x5Top;
   float e2x5Bottom;
   float e2x5Left;
   float e2x5Right;

   bool electronMatch;
};
