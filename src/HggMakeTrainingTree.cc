#include <HggMakeTrainingTree.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"
using namespace std;
#define PI atan(1.0)*4

void HggMakeTrainingTree::processConfig(ReadConfig &cfg){
}

void HggMakeTrainingTree::clear(){
}

void HggMakeTrainingTree::firstInit(){
  setProcessCollection("Electrons",requireGenMatchElectron);
  setProcessCollection("Jets",false);
  setProcessCollection("Muons",false);
}

void HggMakeTrainingTree::processEntry(Long64_t iEntry){  
  if(nVtx==0) return;
  evtNumberOut = evtNumber;
  nVtxOut = nVtx;
  rhoOut = rho;

  for(auto& photonIt : *Photons_) {
    outPhoton = &photonIt;
    
    if(outPhoton->SC.energy/cosh(outPhoton->SC.eta) < minPt) continue;
    if(outPhoton->genMatch.index == -1 && requireGenMatchPhoton) continue;

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

void HggMakeTrainingTree::setupOutputTree(){
    outTree = new TTree("hggPhotonTree","");
    outTree->Branch("evtNumber",&evtNumberOut,"evtNumber/L");
    outTree->Branch("Photon",&outPhoton);
    outTree->Branch("nVtx",&nVtx,"nVtx/I");
    outTree->Branch("rho",&rho,"rho/F");
}
