#include <HggMakeVariableTree.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"
using namespace std;
#define PI atan(1.0)*4

void HggMakeVariableTree::processConfig(ReadConfig &cfg){
  photonID = new HggPhotonID;
  photonID->setConfig(cfg.getPath());
  photonID->Init();
}

void HggMakeVariableTree::clear(){
}

void HggMakeVariableTree::firstInit(){
  setProcessCollection("Jets",false);
  setProcessCollection("Muons",false);
}

void HggMakeVariableTree::processEntry(Long64_t iEntry){  
  if(nVtx==0) return;
  evtNumberOut = evtNumber;
  nVtxOut = nVtx;
  rhoOut = rho;

  photonID->setVertices(nVtx,vtxX,vtxY,vtxZ);
  TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);


  std::array<float,2> pass_pt={-1,-1};
  std::array<TLorentzVector,2> pass_4vec;
  for(auto& photonIt : *Photons_) {
    if( !photonID->getPreSelection(photonIt, nVtx_, rho, 0) ) continue;
    TLorentzVector p4 = photonIt->p4FromVtx(vtx,photonIt->energy,false);
    flaot pt = p4.Pt();
    for(int i=0; i< pass_pt.size(); i++) {
      if(pt > pass_pt[i]) {
	std::swap(pt,pass_pt[i]);
	std::swap(
      }
    }
    
    
    //match Electrons
    if(requireGenMatchElectron) {
      bool match=false;
      for(auto& electronIt: *Electrons_) {
	if(electronIt.SC.index == photonIt.SC.index) {
	  match=true;
	  outPhoton->genMatch = electronIt.genMatch;
	  break;	  
	}
      }

      if(match==false) continue;
    }

    outTree->Fill();
      
  }
  
}

void HggMakeVariableTree::setupOutputTree(){
  outTree = new TTree("output","");

   outTree->Branch("se",&se);
   outTree->Branch("etaSC",&etaSC);
   outTree->Branch("r9",&r9);
   outTree->Branch("phi",&phi);
   outTree->Branch("pt",&pt);
   outTree->Branch("realPho",&realPho);
   outTree->Branch("realEle",&realEle);
   outTree->Branch("passPre",&passPre);
   outTree->Branch("pfChargedGood",&pfChargedGood);
   outTree->Branch("pfChargedWorst",&pfChargedWorst);
   outTree->Branch("pfPhoton",&pfPhoton);
   outTree->Branch("HE",&HE);
   outTree->Branch("sieie",&sieie);
   outTree->Branch("sieip",&sieip);
   outTree->Branch("sipip",&sipip);
   outTree->Branch("passCiC",&passCiC);
   outTree->Branch("passCiC_id",&passCiC_id);
   outTree->Branch("passCiC_iso",&passCiC_iso);
   outTree->Branch("Trigger",&Trigger);
   outTree->Branch("TightPt",&TightPt);
   outTree->Branch("mass",&mass);
   outTree->Branch("e3x3",&e3x3);
   outTree->Branch("e5x5",&e5x5);
   outTree->Branch("rawE",&rawE);
   outTree->Branch("etaWidth",&etaWidth);
   outTree->Branch("phiWidth",&phiWidth);
   outTree->Branch("nBC",&nBC);
   outTree->Branch("energyBC",&energyBC);
   outTree->Branch("etaBC",&etaBC);
   outTree->Branch("phiBC",&phiBC);

   outTree->Branch("eMax",&eMax);
   outTree->Branch("e2nd",&e2nd);
   outTree->Branch("eTop",&eTop);
   outTree->Branch("eBottom",&eBottom);
   outTree->Branch("eLeft",&eLeft);
   outTree->Branch("eRight",&eRight);

   outTree->Branch("e2x5Max",&e2x5Max);
   outTree->Branch("e2x5Top",&e2x5Top);
   outTree->Branch("e2x5Bottom",&e2x5Bottom);
   outTree->Branch("e2x5Left",&e2x5Left);
   outTree->Branch("e2x5Right",&e2x5Right);

   outTree->Branch("electronMatch",&electronMatch);
   outTree->Branch("position",&pos);

}
