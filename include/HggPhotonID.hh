#include <VecbosEGObject.hh>
#include "ReadConfig.hh"

#include <vector>
#include <iostream>
#include <string>
#include <memory>

#include <TChain.h>
#include <TH1F.h>
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TVector3.h"
using namespace std;

class VertexOutOfRange : public std::exception {
public:
  virtual const char * what() const throw() {
    std::cout << "FATAL ERROR: Trying to access a vertex out of range" << std::endl;
  }
};
                                    
class HggPhotonID{
public:
  HggPhotonID();
  ~HggPhotonID();
  void setConfig(std::string s){ configFile = s; }

  bool isValid(){return valid;}
  void Init();

  bool getPreSelection(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex);
  float getIdMVA(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex);
  bool getIdCiCPF(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex);
  bool getIdCiC(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex);
  
  int getCiCCat(VecbosPho* pho);

  bool getPreSelection2011(VecbosPho* pho,int nVertex, float rhoFastJet, int selVtxIndex);
  bool getPreSelectionMay2012(VecbosPho* pho,int nVertex, float rhoFastJet, int selVtxIndex);

  std::map<std::string,TH1F*>* getHists(){return &InputDists;}

  void setVertices(int nPV, float* xPV, float* yPV, float* zPV);

  void fillIsoVariables(VecbosPho* pho, ReducedPhotonData* data, int nVertex, float rhoFastJet, int selVtxIndex);
  
  void setDoEcalIso(bool b){doECALIso=b;}
private:
  std::string configFile;
  bool valid;

  bool doECALIso;

  string version;

  static constexpr float isoSumConst = 0;
  static constexpr float isoSumConstPF = 2.5;
  static constexpr float rhoFac = 0.09;
  static constexpr float rhoFacBad = 0.23;

  string weightFile_IdEB_2011;
  string weightFile_IdEE_2011;
  string weightFile_IdEB_2012;
  string weightFile_IdEE_2012;
  string methodName_Id;

  TMVA::Reader *photonMVA_EB_2011;
  TMVA::Reader *photonMVA_EE_2011;
  TMVA::Reader *photonMVA_EB_2012;
  TMVA::Reader *photonMVA_EE_2012;
  void setupTMVA();

  std::map<std::string,TH1F*> InputDists;

  void fillVariables(VecbosPho* pho,int nVertex, float rhoFastJet, int selVtxIndex);
  //photon ID MVA
  float hoe;
  float sigietaieta;
  float isosumoet;
  float isosum;
  float trkisooet;
  float isosumoetbad;
  float isosumbad;
  float r9;
  float ecalisodr03;
  float ecalisodr04;
  float hcalisodr03;
  float hcalisodr04;
  float nVertexf;
  float etasc;
  float scetawidth;
  float scphiwidth;

  float pfChargedIsoGood03;
  float pfChargedIsoBad04;
  float pfChargedIsoGood03oet;
  float pfChargedIsoBad04oet;
  float pfChargedIsoZero03;
  float pfChargedIsoZero03oet;
  float pfPhotonIso03;
  float pfPhotonIso03oet;
  float pfPhotonIso04;
  float pfPhotonIso04oet;
  float sigietaiphi;
  float s4Ratio;
  float rho;
  float sigRR;
  float isosumoetPF;
  float isosumoetbadPF;
  float isosumPF;
  float isosumbadPF;

  float eT;
  float eTBad;
  std::vector<TVector3> vertices;  
  TVector3 selVtxPos;
};
