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

void MuMuGammaSelector::processConfig(ReadConfig &cfg){
  photonID = new HggPhotonID;
  photonID->setConfig(cfg.getPath());
  photonID->Init();
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

void MuMuGammaSelector::firstInit(){
  //setProcessCollection("Electrons",false);
  //setProcessCollection("Jets",false);
}

void MuMuGammaSelector::processEntry(Long64_t iEntry){  
  if(nVtx==0) return;
  this->writeEventInfo();
  
  photonID->setVertices(nVtx,vtxX,vtxY,vtxZ);
  
  TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);
  
  for(int iMu=0;iMu<nMu_;iMu++){
    VecbosMu mu1 = Muons_->at(iMu);
    if(!passPresel(mu1)) return;
    for(int jMu=iMu+1; jMu<nMu_;jMu++){
      VecbosMu mu2 = Muons_->at(jMu);
      if(!passPresel(mu2)) return;
      
      if(mu1.pt < 20 && mu2.pt < 20) return;
      if(mu1.charge*mu2.charge >=0) return;
      
      for(int iPho=0;iPho<nPho_;iPho++){
	VecbosPho pho = Photons_->at(iPho);
	
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
}
void MuMuGammaSelector::writeEventInfo(){
  runNumberOut = runNumber;
  evtNumberOut = evtNumber;
  nVtxOut = nVtx;
  rhoOut = rho;
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






