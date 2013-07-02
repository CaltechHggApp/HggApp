#include <MuMuGammaSelector.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"
using namespace std;
#define PI atan(1.0)*4

MuMuGammaSelector::MuMuGammaSelector(vector<string> fNames, string treeName,string outFName):
  fChain(0),
  isData_(true),
  __Muons(0),
  __Photons(0)
{
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

MuMuGammaSelector::~MuMuGammaSelector(){
  delete fChain;
}


void MuMuGammaSelector::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int MuMuGammaSelector::init(){
  if(!valid) return -1;
  
  this->setBranchAddresses();
  this->setupOutputTree();

  photonID = new HggPhotonID;
  photonID->setConfig(cfg);
  photonID->Init();
  return 0;
}

void swap(VecbosEle& e1, VecbosEle &e2){
  VecbosEle t= e1; e1=e2; e2=e1;
}

void MuMuGammaSelector::clear(){
  nMuMuG=0;
  outMuon1.clear();
  outMuon2.clear();
  outPhoton.clear();
}

void MuMuGammaSelector::Loop(){
  if(!valid) return;

  this->init();
  cout << "Getting Entries ... "  << endl;
  int nEntries = fChain->GetEntries();
  int nbytes = 0;

  Long64_t iEntry=-1;
  while(fChain->GetEntry(++iEntry)){
    if(nVtx==0) continue;
    this->clear();
    this->writeEventInfo();

    photonID->setVertices(nVtx,vtxX,vtxY,vtxZ);

    TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);

    for(int iMu=0;iMu<__nMu;iMu++){
      VecbosMu mu1 = __Muons->at(iMu);
      if(!passPresel(mu1)) continue;
      for(int jMu=iMu+1; jMu<__nMu;jMu++){
	VecbosMu mu2 = __Muons->at(jMu);
	if(!passPresel(mu2)) continue;
	
	if(mu1.pt < 20 && mu2.pt < 20) continue;
	if(mu1.charge*mu2.charge >=0) continue;

	for(int iPho=0;iPho<__nPho;iPho++){
	  VecbosPho pho = __Photons->at(iPho);
	  
	  //fill the different masses
	  TLorentzVector p4Mu1; p4Mu1.SetPtEtaPhiM(mu1.pt,mu1.eta,mu1.phi,0.106);
	  TLorentzVector p4Mu2; p4Mu2.SetPtEtaPhiM(mu2.pt,mu2.eta,mu2.phi,0.106);

	  massMuMu[nMuMuG] = (p4Mu1+p4Mu2).M();
	  massMuMuGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.energy,false)).M();
	  massMuMuRegGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.correctedEnergy,false)).M();
	  massMuMuScaleGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.scaledEnergy,false)).M();
	  if(pho.genMatch.index>=0) massMuMuGenGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.genMatch.energy,false)).M();
	  else massMuMuGenGamma[nMuMuG] = -1;
	  
	  outMuon1.push_back(mu1);
	  outMuon2.push_back(mu2);
	  outPhoton.push_back(pho);
	  //float eT = pho.p4FromVtx(vtx,pho.energy,false).Et();
	  isosumoetPho[nMuMuG] = 0;
	  mvaPho[nMuMuG] = photonID->getIdMVA(&pho,nVtx,rho,0); 
	
	  nMuMuG++;

	}//for(int iPho=0;
	
      }//for(int jMu=iMu+1;

    }//for(int iMu=0;
    outTree->Fill();

  }
  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outTree->Write();
  f->Close();
}

void MuMuGammaSelector::writeEventInfo(){
  runNumberOut = runNumber;
  evtNumberOut = evtNumber;
  nVtxOut = nVtx;
  rhoOut = rho;
}


void MuMuGammaSelector::setBranchAddresses(){
  if(!valid) return;
  //Event info
  fChain->SetBranchAddress("runNumber",&runNumber);
  fChain->SetBranchAddress("evtNumber",&evtNumber);
  //fChain->SetBranchAddress("isRealData",&_isData);
  
  //information for the Electrons
  fChain->SetBranchAddress("nMu",&__nMu);
  fChain->SetBranchAddress("Muons",&__Muons);
  fChain->SetBranchAddress("nPho",&__nPho);
  fChain->SetBranchAddress("Photons",&__Photons);
 
  // information for the vertex
  fChain->SetBranchAddress("nVtx", &nVtx);
  fChain->SetBranchAddress("vtxX",vtxX);
  fChain->SetBranchAddress("vtxY",vtxY);
  fChain->SetBranchAddress("vtxZ",vtxZ);
  fChain->SetBranchAddress("rho", &rho);

}

void MuMuGammaSelector::setupOutputTree(){
    outTree = new TTree("MuMuGamma","");
    outTree->Branch("runNumber",&runNumberOut);
    outTree->Branch("evtNumber",&evtNumberOut);
    outTree->Branch("nVtx",&nVtxOut);
    outTree->Branch("rho",&rhoOut);

    outTree->Branch("nMuMuG",&nMuMuG);
    outTree->Branch("massMuMuGamma",massMuMuGamma,"massMuMuGamma[nMuMuG]");
    outTree->Branch("massMuMuRegGamma",massMuMuRegGamma,"massMuMuRegGamma[nMuMuG]");
    outTree->Branch("massMuMuScaleGamma",massMuMuScaleGamma,"massMuMuScaleGamma[nMuMuG]");
    outTree->Branch("massMuMuGenGamma",massMuMuGenGamma,"massMuMuGenGamma[nMuMuG]");
    outTree->Branch("massMuMu",massMuMu,"massMuMu[nMuMuG]");
    outTree->Branch("puWeight",puWeight,"puWeight[nMuMuG]");
    outTree->Branch("Muon1",&outMuon1);
    outTree->Branch("Muon2",&outMuon2);
    outTree->Branch("Photon",&outPhoton);
    outTree->Branch("isosumoetPho",isosumoetPho,"isosumoetPho[nMuMuG]");
    outTree->Branch("mvaPho",mvaPho,"mvaPho[nMuMuG]");
}

bool MuMuGammaSelector::passPresel(VecbosMu &mu){
  if(mu.pt < 10) return false;
  if(!mu.isGlobalMuon || !mu.isTrackerMuon) return false;
  if(mu.nTrackHits <= 10 || mu.nPixelHits==0) return false;
  if(mu.trackImpactPar >=0.2) return false;
  if(mu.trkIso >= 3) return false;
 
  return true;
}






