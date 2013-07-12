#include <BaseSelector.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TRandom3.h"
using namespace std;

#define debugSelector 0

BaseSelector::BaseSelector():
  fChain(0),
  valid(false),
  isData_(true),
  outFile(0),
  outTree(0),
  Photons_(0),
  ggVerticesPhotonIndices(0),
  ggVerticesVertexIndex(0),
  ggVerticesPerEvtMVA(0),
  Muons_(0),
  Electrons_(0),
  Jets_(0),
  GenHiggs(0),
  GenPhotons(0)
{
}

BaseSelector::BaseSelector(vector<string> fNames, string treeName,string outFName):
  fChain(0),
  valid(false),
  isData_(true),
  outFile(0),
  outTree(0),
  Photons_(0),
  ggVerticesPhotonIndices(0),
  ggVerticesVertexIndex(0),
  ggVerticesPerEvtMVA(0),
  Muons_(0),
  Electrons_(0),
  Jets_(0),
  GenHiggs(0),
  GenPhotons(0)
{
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

void BaseSelector::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int BaseSelector::baseInit(){
  if(!valid) return -1;

  triggerDec = new int[triggers.size()];

  
  this->setBranchAddresses();
  this->setupOutputTree();

  ReadConfig cfg;
  int retcode = cfg.read(configFile);
  if(retcode != 0){
    cout << "Error reading configuration file!" <<std::endl
	 << "Error Code: " << retcode << std::endl;
    valid = false;
    return -1;
  }
  processConfig(cfg);

  outFile = new TFile(outputFile.c_str(),"RECREATE");

  int slaveReturn = init(); // run the initialization of the derived class
  if(slaveReturn!=0){
    outFile->Close();
  }
  return slaveReturn;
}

void BaseSelector::Loop(){
  if(!valid) return;


  //run ths initialization
  if(this->baseInit()!=0){
    std::cout << "ERROR INITIALIZING ... ABORTING!" <<std::endl;
    return;
  }
  this->setDefaults();

  cout << "Getting Entries ... "  << endl;
  Long64_t nEntries = fChain->GetEntries();
  Long64_t jentry=-1;

  while(fChain->GetEntry(++jentry)){
    if(jentry%500==0) cout << ">> Processing Entry " << jentry << "/" << nEntries << endl;


    clear();
    setDefaults();
    processEntry(jentry);
    
    outTree->Fill();
  }//while(fChain...

  outFile->cd();
  outTree->Write();
  write();
  outFile->Close();
}

void BaseSelector::setBranchAddresses(){
  if(!valid) return;
  //Event info
  fChain->SetBranchAddress("lumiBlock",&lumiBlock);
  fChain->SetBranchAddress("runNumber",&runNumber);
  fChain->SetBranchAddress("evtNumber",&evtNumber);
  //fChain->SetBranchAddress("isRealData",&_isData);
  
  fChain->SetBranchAddress("eeBadScFilterFlag",&eeBadScFilterFlag);
  fChain->SetBranchAddress("hcalLaserEventFilterFlag",&hcalLaserEventFilterFlag);
  fChain->SetBranchAddress("HBHENoiseFilterResultFlag",&HBHENoiseFilterResultFlag);
  fChain->SetBranchAddress("isNotDeadEcalCluster",&isNotDeadEcalCluster);
  fChain->SetBranchAddress("trackerFailureFilterFlag",&trackerFailureFilterFlag);
  fChain->SetBranchAddress("CSCHaloFilterFlag",&CSCHaloFilterFlag);
  fChain->SetBranchAddress("drDead",&drDead);
  fChain->SetBranchAddress("drBoundary",&drBoundary);
  fChain->SetBranchAddress("ECALTPFilterFlag",&ECALTPFilterFlag);

 ///information for the vertex
  fChain->SetBranchAddress("nVtx",&nVtx);
  fChain->SetBranchAddress("vtxX",vtxX);
  fChain->SetBranchAddress("vtxY",vtxY);
  fChain->SetBranchAddress("vtxZ",vtxZ);
  fChain->SetBranchAddress("vtxChi2",vtxChi2);
  fChain->SetBranchAddress("vtxNdof",vtxNdof);
  fChain->SetBranchAddress("vtxNormalizedChi2",vtxNormalizedChi2);
  fChain->SetBranchAddress("vtxTrackSize",vtxTrackSize);
  fChain->SetBranchAddress("vtxIsFake",vtxIsFake);
  fChain->SetBranchAddress("vtxIsValid",vtxIsValid);
  
  fChain->SetBranchAddress("rho", &rho);
  fChain->SetBranchAddress("rhoEtaMax44", &rhoEtaMax44);

  fChain->SetBranchAddress("pileupWeight", &pileupWeight);
 
  //objects
  fChain->SetBranchAddress("nPho",&nPho_);
  fChain->SetBranchAddress("Photons",&Photons_);
  
  fChain->SetBranchAddress("photonMatchedElectron",photonMatchedElectron);
  fChain->SetBranchAddress("nPair",&nPair_); 
  fChain->SetBranchAddress("ggVerticesPhotonIndices",&ggVerticesPhotonIndices);
  fChain->SetBranchAddress("ggVerticesVertexIndex",&ggVerticesVertexIndex);
  fChain->SetBranchAddress("ggVerticesPerEvtMVA",&ggVerticesPerEvtMVA);

  fChain->SetBranchAddress("nMu",&nMu_);
  fChain->SetBranchAddress("Muons",&Muons_);
  
  fChain->SetBranchAddress("nJet",&nJet_);
  fChain->SetBranchAddress("Jets",&Jets_);

  fChain->SetBranchAddress("nGenHiggs",&nGenHiggs);
  fChain->SetBranchAddress("GenHiggs",&GenHiggs);
  
  fChain->SetBranchAddress("nGenPho",&nGenPho);
  fChain->SetBranchAddress("GenPhotons",&GenPhotons);
  
  fChain->SetBranchAddress("nPU",&inPU);
  
  fChain->SetBranchAddress("PFMET",&pfMet);
  fChain->SetBranchAddress("PFMETPhi",&pfMetPhi);

  vector<string>::const_iterator trigIt;
  int i=0;
  for(trigIt=triggers.begin();trigIt!=triggers.end();trigIt++,i++){
    fChain->SetBranchAddress(trigIt->c_str(),&(triggerDec[i]));
  }

}

bool BaseSelector::requireTrigger(){
  if(!_isData) return true; //no triggers on MC

  if(triggers.size()==0) return true; //no trigger selection
  
  for(int i=0;i<triggers.size();i++){
    if(triggerDec[i]) return true;
  }
  return false;
}
