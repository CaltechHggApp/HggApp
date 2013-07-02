 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef HggVertexing_h
#define HggVertexing_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<unsigned short> >+; 
#endif
//#include "TMVAGui.C"                         
//#if not defined(__CINT__) || defined(__MAKECINT__) 
//#include "../h2ginclude/HggVertexAnalyzer.h"
//#include "../h2ginclude/HggVertexFromConversions.h"
//#include "../h2ginclude/PhotonInfo.h"
//#include "../h2ginclude/VertexAlgoParameters.h"
                 
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TVector3.h"
//#endif 
using namespace TMVA; 
using namespace std;
#include "VecbosEGObject.hh"


class HggVertexing{
public:
  HggVertexing(VecbosBase*); /// Class Constructor
  virtual ~HggVertexing();     /// Class Destructor
  void setConfigFile(string s){configFilePath = s;}

  std::vector<pair<int,float> > vertex_tmva(VecbosPho*,VecbosPho*,float&);
  void useConversions(bool b=true){useConversion = b;}
  void init(); // do variable initialization 
  bool getIsInit(){return isInit;}
  void saveInputs(TTree*);
  void doRescaleTrkPt(bool b){rescaleTrkPt = b;}

private:
  VecbosBase *base;
  bool isInit;
  bool doSaveInputs;
  bool rescaleTrkPt;

  std::pair<float,float> getZConv(VecbosPho*,VecbosPho*);
  std::pair<float,float> getZConv(VecbosPho*);
  float convCorrectedDz(VecbosConversion*,TVector3);
  float Z0EcalVtxCiC(VecbosConversion*,TVector3,TVector3);
  int getNConv(VecbosPho* p1, VecbosPho* p2);
  bool isGoodVertex(int);

  TMVA::Reader * perVtxReader; 
  TMVA::Reader * perEvtReader;
  float *PerVtxVars;
  float *PerEvtVars; //variables for the tmva


  std::vector<std::pair<int,float> > evalPerVtxMVA(VecbosPho*, VecbosPho*);
  float evalPerVtxMVA(float ptbal, float ptasym, float logsumpt2, float limPullToConv, float nConv);
  float evalPerEvtMVA(VecbosPho* pho1,VecbosPho* pho2,std::vector<std::pair<int,float> > *perVertexRank);

  //read these from the config file
  string configFilePath;
  string perVtxMvaWeights;
  string perVtxMvaMethod;
  string perEvtMvaWeights; 
  string perEvtMvaMethod;

  //options
  bool useConversion;
};
#endif


