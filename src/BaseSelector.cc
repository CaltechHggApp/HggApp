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
  Jets_(0)
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
  Jets_(0)
{
  std::cout << outFName << std::endl;
  outFile = new TFile(outFName.c_str(),"RECREATE");
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

void BaseSelector::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    std::cout << "chaining: " << *name << std::endl;
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int BaseSelector::baseInit(){
  if(!valid) return -1;

  firstInit();

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
    
    if(doFill) outTree->Fill();
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
 
  fChain->SetBranchStatus("*",1);
  for(auto& colIt: CollectionsToProcess){
    if(colIt.second==false){
      fChain->SetBranchStatus(colIt.first,0);
    }
  }

  //objects
  if(CollectionsToProcess["Photons"]){
    fChain->SetBranchAddress("nPho",&nPho_);
    fChain->SetBranchAddress("Photons",&Photons_);
    
    fChain->SetBranchAddress("photonMatchedElectron",photonMatchedElectron);
    fChain->SetBranchAddress("nPair",&nPair_); 
    fChain->SetBranchAddress("ggVerticesPhotonIndices",&ggVerticesPhotonIndices);
    fChain->SetBranchAddress("ggVerticesVertexIndex",&ggVerticesVertexIndex);
    fChain->SetBranchAddress("ggVerticesPerEvtMVA",&ggVerticesPerEvtMVA);
  }
  if(CollectionsToProcess["Muons"]){
    fChain->SetBranchAddress("nMu",&nMu_);
    fChain->SetBranchAddress("Muons",&Muons_);
  }
  if(CollectionsToProcess["Electrons"]){
    fChain->SetBranchAddress("nEle",&nEle_);
    fChain->SetBranchAddress("Electrons",&Electrons_);
  }
  if(CollectionsToProcess["Jets"]){
    fChain->SetBranchAddress("nJet",&nJet_);
    fChain->SetBranchAddress("Jets",&Jets_);
  }
  
  fChain->SetBranchAddress("nGenHiggs",&nGenHiggs);
  fChain->SetBranchAddress("GenHiggs",&GenHiggs);
  
  fChain->SetBranchAddress("nGenPho",&nGenPho);
  fChain->SetBranchAddress("GenPhotons",&GenPhotons);
  
  fChain->SetBranchAddress("nGenEle",&nGenEle);
  fChain->SetBranchAddress("GenElectrons",&GenElectrons);
  
  fChain->SetBranchAddress("nGenMu",&nGenMu);
  fChain->SetBranchAddress("GenMuons",&GenMuons);
  
  fChain->SetBranchAddress("nGenOthers",&nGenOthers);
  fChain->SetBranchAddress("GenOthers",&GenOthers);
  
  fChain->SetBranchAddress("nPU",&inPU);
  
  fChain->SetBranchAddress("PFMET",&pfMet);
  fChain->SetBranchAddress("PFMETPhi",&pfMetPhi);

  fChain->SetBranchAddress("type1PFMET",&type1PfMet);
  fChain->SetBranchAddress("type1PFMETPhi",&type1PfMetPhi);

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
