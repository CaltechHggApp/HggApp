#include <VecbosEGObject.hh>
#include <vector>
#include <iostream>
#include <exception>
#include <string>

#include <TChain.h>

#include "ReadConfig.hh"
using namespace std;


class BaseSelector{
public:
  BaseSelector();
  BaseSelector(vector<string> fNames,string treeName,string outputFile);
  void loadChain(vector<string> fNames, string treeName);
  void setOutputFile(string s){outputFile = s;}
  void setConfig(string s){configFile = s;}
  bool isValid(){return valid;}
  void addTrigger(string s){triggers.push_back(s);}
  void setIsData(bool d){isData_=d;}
  void Loop();
  void setDoFill(bool b){doFill=b;}
protected:
  bool doFill=true; // should the main event loop do the filling (true)


  TChain* fChain;
  bool valid;
  bool isData_;

  TFile *outFile;
  TTree* outTree;

  string configFile;
  string outputFile;
  
  std::vector<std::string> triggers;
  int *triggerDec;
  bool requireTrigger();

  //implemneted by the base class only
  int baseInit();
  void setBranchAddresses();



  //the following must be implemented
  virtual int init()=0;
  virtual void processEntry(Long64_t iEntry)=0;
  virtual void setupOutputTree()=0;

  //the following may be implemented
  virtual void processConfig(ReadConfig& config){}
  virtual void setDefaults(){}
  virtual void clear(){}
  virtual void write(){}
  virtual void firstInit(){}

  //We can turn on and off different inputs here:
  std::map<TString,bool> CollectionsToProcess = {
    { "Photons",true },
    { "Jets", true},
    { "Muons", true },
    { "Electrons", true}
  };
  void setProcessCollection(TString n, bool b){ 
    try{
      CollectionsToProcess.at(n) = b;
    }catch(std::exception& e){
      std::cout << "FATAL ERROR:  Trying to set the status of invalid collection " << n << std::endl;
      std::cout << "Valid Collections: " << std::endl;
      for(auto& it : CollectionsToProcess) {
	std::cout << it.first << "    process: " << it.second << std::endl;
      }
      throw e;
    }
  }
     

  //
  // INPUT VARIABLES
  //
  const static int maxPho = 100;
  int nPho_;
  std::vector<VecbosPho> *Photons_; // this contains ALL photons
  
  bool photonMatchedElectron[maxPho];

  //this is the collection of the TMVA selected vertices for each photon pair
  //the format is std::pair< (iPho1 << 14) + iPho2,iVrt> 
  int nPair_;
  std::vector<std::pair<int,int> >     *ggVerticesPhotonIndices;
  std::vector<std::pair<int, float> >  *ggVerticesVertexIndex;
  std::vector<float>                   *ggVerticesPerEvtMVA;

  int nMu_;
  std::vector<VecbosMu> *Muons_;

  int nEle_;
  std::vector<VecbosEle> *Electrons_;

  int nJet_;
  std::vector<VecbosJet> *Jets_;
  
  // for each photon, a vector of floats giving the track iso from each ggVertex

  int nVtx; 
  static const int MAXVX = 100;
  float vtxX[MAXVX];
  float vtxY[MAXVX];
  float vtxZ[MAXVX];
  float vtxChi2[MAXVX];
  float vtxNdof[MAXVX];
  float vtxNormalizedChi2[MAXVX];
  int vtxTrackSize[MAXVX];
  int vtxIsFake[MAXVX];
  int vtxIsValid[MAXVX];

  float rho;
  float rhoEtaMax44;

  int lumiBlock;
  int runNumber;
  long evtNumber;
  bool _isData;

  float evtWeight;
  float pileupWeight;
  
  int nGenHiggs;
  std::vector<VecbosGen> *GenHiggs;

  int nGenPho;
  std::vector<VecbosGen> *GenPhotons;

  float pfMet;
  float pfMetPhi;

  float inPU;

  bool eeBadScFilterFlagOut;
  bool hcalLaserEventFilterFlagOut;
  bool HBHENoiseFilterResultFlagOut;
  bool isNotDeadEcalClusterOut;
  bool trackerFailureFilterFlagOut;
  bool CSCHaloFilterFlagOut;
  bool drDeadOut; 
  bool drBoundaryOut;
  bool ECALTPFilterFlagOut;

  bool eeBadScFilterFlag;
  bool hcalLaserEventFilterFlag;
  bool HBHENoiseFilterResultFlag;
  bool isNotDeadEcalCluster;
  bool trackerFailureFilterFlag;
  bool CSCHaloFilterFlag;
  bool drDead; 
  bool drBoundary;
  bool ECALTPFilterFlag;
};
