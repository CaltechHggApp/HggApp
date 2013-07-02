#ifndef HggMakePhotonTree_hh
#define HggMakePhotonTree_hh

#include <VecbosBase.hh>
#include <TTree.h>
#include <TChain.h>
#include <string>
#include <TLorentzVector.h>

struct ECAL_GEO;

class HggMakePhotonTree : public VecbosBase{
public:
  HggMakePhotonTree(TTree*,std::string,bool);
  void Loop(int,int);
private:
  const static float rhoFac    = 0.17;
  const static float rhoFacBad = 0.52;
  const static float isoSumConst = 0;//5;
  const static float isoSumConstBad = 0;//7;
  bool usePF;
  string outFName;
  void init();
  void setupBranches();
  float computeTrackIso(int,int,float,float,float,float,float,float);
  int genMatchPho();
  int convMatchPho();
  TLorentzVector p4FromVtx(int,int,float);
  void fillGenInfo();
  void fillGapVariables();
  void fillIsoVariables();
  void fillConvVariables();
  void fillEnergyVariables();

  ECAL_GEO ecalGeometry;

  int PhoIndex;
  int GenPhoIndex;
  int ConvPhoIndex;
  int SCPhoIndex;
  float phoPt;

  TTree *tree;
  //variables for the tree:
  //keep the naming the same as Yong's code ....
  float etaCGap;
  float phiCGap;
  float etaSGap;
  float phiSGap;
  float etaMGap;
  float phiMGap;

  float isosumoet;
  float isosumoetbad;
  float trkisooet;
  float hoe;
  float r9;
  int haspromptele;
  float ecalisodr03;
  float ecalisodr04;
  float hcalisodr03;
  float hcalisodr04;
  float trkisodr03vtxbestGenMatched;
  float trkisodr04vtxWorst;

  float etrue;
  float etatrue;
  float phitrue;
  float etaecaltrue;
  float phiecaltrue;

  float escraw;
  float e5x5;
  float eps;
  float esc;
  float epht;
  float etasc;
  float phisc;

  float e2x2;
  float e3x3;
  float e1x3;
  float e3x1;
  float e4x4;
  float e2x5;
  float e2x5right;
  float e2x5left;
  float e2x5top;
  float e2x5bottom;
  float e2x5max;

  float etapht;
  float phipht;
  
  float scetawidth;
  float scphiwidth;
  float sigietaieta;
  float sigiphiiphi;
  float sigcovietaiphi;
  
  float emax;
  int scnbc;
  int scncrystal;
  float scbcfmax;
  int seedieta;
  int seediphi;
  float eleft;
  float eright;
  float etop;
  float ebottom;
  float e2max;

  float convp;
  float convz;
  float convrho;
  float convpttrk1;
  float convpttrk2;
  float convpt;
  float deltaphi_convsc;
  float deltaeta_convsc;
  int truephtmatched;

  int pidphtmom;

  float rho;
  int nVertex;

  int lumiBlock;
  int runNumber;
  int evtNumber;
  int isRealData;

  float dbc2bceta;
  float dbc2bcphi;

  float bce2nd;
  float dbceta;
  float dbcphi;
  float bce;
  float bc2e;
  float bc2emax;
  float bc2e2nd;
  float bc2eleft;
  float bc2eright;
  float bc2etop;
  float bc2ebottom;
  float dbc2eta;
  float dbc2phi;
  float bc2sigietaieta;
  float bc2sigiphiiphi;
  float bc2sigietaiphi;
  float bc2e3x3;
  float bc2e5x5;
  
  float bclaste;
  float bclaste3x3;
  float bclaste5x5;
  float bclastsigietaieta;
  float bclastsigiphiiphi;
  float bclastsigietaiphi;
  float dbclasteta;
  float dbclastphi;
  
  float bclast2e;
  float bclast2e3x3;
  float bclast2e5x5;

  float bclast2sigietaieta;
  float bclast2sigiphiiphi;
  float bclast2sigietaiphi;
  
  float dbclast2eta;
  float dbclast2phi;
  float bcseedetacry;
  
  float bc2etacry;
  float bclast2etacry;
  float bclastetacry;
  float bcseedphicry;
  float bc2phicry;
  float bclast2phicry;
  float bclastphicry;
  int bc2ieta;
  int bc2iphi;
};


#endif
