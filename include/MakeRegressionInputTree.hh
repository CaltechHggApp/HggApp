#ifndef MakeRegressionInputTree_hh
#define MakeRegressionInputTree_hh

#include "BaseSelector.hh"
#include "StandardPhotonID.hh"

class MakeRegressionInputTree: public BaseSelector {
public:
  MakeRegressionInputTree(std::vector<std::string> fNames, std::string treeName,std::string outputFile);
  
private:

  StandardPhotonID photonID;

  //mandatory overrides
  virtual void processEntry(Long64_t iEntry) override;
  virtual int init() override;
  virtual void processConfig(ReadConfig& cfg) override {}
  virtual void setupOutputTree() override;

  virtual void clear() override;

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
  int pu;

  int numVtx;
  float jetRho;

  float pfChargedGood;
  float pfChargedWorst;
  float pfPhoton;
  float HE;
  float sieie;
  float sieip;
  float sipip;
  
  int passWP90;
  int passWP90_id;
  int passWP90_iso;
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

#endif
