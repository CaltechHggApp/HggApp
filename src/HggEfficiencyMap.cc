#include <HggEfficiencyMap.hh>
#include "ReadConfig.hh"
#include "CommonTools/include/Utils.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"
using namespace std;
using namespace TMVA;

#define debugSelector 0

HggEfficiencyMap::HggEfficiencyMap():
  fChain(0),
  valid(false),
  Photons_(0),
  Muons_(0),
  Jets_(0),
  GenHiggs(0),
  GenPhotons(0),
  GenElectrons(0),
  mode(kMC)
{
}

HggEfficiencyMap::HggEfficiencyMap(vector<string> fNames, string treeName,string outFName):
  fChain(0),
  valid(false),
  Photons_(0),
  Muons_(0),
  Jets_(0),
  GenHiggs(0),
  GenPhotons(0),
  GenElectrons(0),
  mode(kMC)
{
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

HggEfficiencyMap::~HggEfficiencyMap(){
  delete fChain;
}


void HggEfficiencyMap::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int HggEfficiencyMap::init(){
  if(!valid) return -1;

  outputTree = new TTree("EfficiencyTree","");
  this->setOutputBranches();
  
  this->setBranchAddresses();


  leadPhoEtMin = 33.;
  subleadPhoEtMin = 25.;
  PtOverMLead = 40./120.;
  PtOverMSubLead = 30./120.;

  PhotonID = new HggPhotonID();
  PhotonID->setConfig(configFile);
  PhotonID->Init();
  return 0;
}

void HggEfficiencyMap::setOutputBranches(){
  outputTree->Branch("GenPt",&outPt);
  outputTree->Branch("RecoPt",&outRecoPt);
  outputTree->Branch("GenEta",&outEta);
  outputTree->Branch("RecoEta",&outRecoEta);
  outputTree->Branch("GenEtaSC",&outEtaSC);
  outputTree->Branch("RecoEtaSC",&outRecoEtaSC);
  outputTree->Branch("GenPhi",&outPhi);
  outputTree->Branch("GenR9",&outR9);

  outputTree->Branch("passMatch",&passMatch,"passMatch/O");
  outputTree->Branch("passEta",&passEta,"passEta/O");
  outputTree->Branch("passPreselection",&passPreselection,"passPreselection/O");
  outputTree->Branch("passCiC",&passCiC,"passCiC/O");
}

void HggEfficiencyMap::Loop(){
  if(!valid) return;

  if(this->init()!=0){
    std::cout << "ERROR INITIALIZING ... ABORTING!" <<std::endl;
    return;
  }

  cout << "Getting Entries ... "  << endl;
  Long64_t nEntries = fChain->GetEntries();
  Long64_t jentry=-1;
  while(fChain->GetEntry(++jentry)){
    if(jentry%500==0) cout << ">> Processing Entry " << jentry << "/" << nEntries << endl;

    this->clear();
    if(debugSelector) cout << "Setting Photon ID Vertices " << endl;
    PhotonID->setVertices(nVtx,vtxX,vtxY,vtxZ);

    if(nVtx==0) continue;
    
    switch(mode){
    case kMC:
      runGen();
      break;

    case kMMG:
      runMMG();
      break;

    case kZeeMC:
      runZeeMC();
      break;

    case kZee:
      runZee();
      break;
    }
    
  }//while(fChain...

  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outputTree->Write();
  f->Close();
}

void HggEfficiencyMap::runGen(){
  //run the MC-based selection
    //DO
    for(int iGen=0; iGen<nGenPho;iGen++){
      VecbosGen& gen = GenPhotons->at(iGen);
      this->clear();
      if(gen.status!=1) continue;
      outPt = gen.pt;
      outEta = gen.eta;
      outPhi = gen.phi;
      //compute the detector eta
      TVector3 vtxPos(vtxX[0],vtxY[0],vtxZ[0]); //gen.Vx,gen.Vy,gen.Vz); // Something is broken with the gen vtx
      TVector3 physVec;
      physVec.SetPtEtaPhi(outPt,outEta,outPhi);
      TVector3 caloVec = (physVec+vtxPos).Unit()*gen.energy;
      TLorentzVector caloP4;
      caloP4.SetVectM(caloVec,0);
      outEtaSC = caloP4.Eta();

      outRecoPt=-1;
      outRecoEta=-99;
      outRecoEtaSC=-99;
      outR9=0;
      VecbosPho* matched = findGenMatchPhoton(&gen);
      if(matched == 0){
	outputTree->Fill(); continue;
      }
      passMatch=true;
      outRecoEta = matched->eta;
      outRecoPt = matched->p4FromVtx(vtxPos,matched->finalEnergy).Pt();
      outRecoEtaSC = matched->SC.eta;
      outR9 = matched->SC.r9;

      if(!etaSelectPhoton(matched)){
	outputTree->Fill(); continue;
      }
      passEta=true;
      
      if(!PhotonID->getPreSelection(matched,nVtx,rho,0)){
	outputTree->Fill(); continue;
      }
      passPreselection=true;

      if(!PhotonID->getIdCiCPF(matched,nVtx,rho,0)){
	outputTree->Fill(); continue;
      }
      passCiC=true;
      
      outputTree->Fill();
    }

    //fill fakes
    for(int i=0;i<nPho_;i++){
      this->clear();
      VecbosPho& pho = Photons_->at(i);
      if(pho.genMatch.index != -1) continue; //only looking for fakes
      bool matched = false;
      outPt = -1;
      outEta = 0;
      outEtaSC = 0;

      TVector3 vtxPos(vtxX[0],vtxY[0],vtxZ[0]); 
      outRecoPt=pho.p4FromVtx(vtxPos,pho.finalEnergy).Pt();
      outRecoEta=pho.eta;
      outRecoEtaSC=pho.SC.eta;
      outR9=pho.SC.r9;
      outPhi = pho.phi;

      
      if(!etaSelectPhoton(&pho)){
	outputTree->Fill(); continue;
      }
      passEta=true;
      
      if(!PhotonID->getPreSelection(&pho,nVtx,rho,0)){
	outputTree->Fill(); continue;
      }
      passPreselection=true;

      if(!PhotonID->getIdCiCPF(&pho,nVtx,rho,0)){
	outputTree->Fill(); continue;
      }
      passCiC=true;
      
      outputTree->Fill();
      
    }
}

void HggEfficiencyMap::runZeeMC(){
  //run the MC-based selection
    //DO
    for(int iGen=0; iGen<nGenEle;iGen++){
      VecbosGen& gen = GenElectrons->at(iGen);
      this->clear();
      if(gen.status!=1) continue;
      outPt = gen.pt;
      outEta = gen.eta;
      outPhi = gen.phi;
      //compute the detector eta
      TVector3 vtxPos(vtxX[0],vtxY[0],vtxZ[0]); //gen.Vx,gen.Vy,gen.Vz); // Something is broken with the gen vtx
      TVector3 physVec;
      physVec.SetPtEtaPhi(outPt,outEta,outPhi);
      TVector3 caloVec = (physVec+vtxPos).Unit()*gen.energy;
      TLorentzVector caloP4;
      caloP4.SetVectM(caloVec,0);
      outEtaSC = caloP4.Eta();

      outRecoPt=-1;
      outRecoEta=-99;
      outRecoEtaSC=-99;
      outR9=0;
      VecbosPho* matched = findGenMatchElectron(&gen);
      if(matched == 0){
	outputTree->Fill(); continue;
      }
      passMatch=true;
      outRecoEta = matched->eta;
      outRecoPt = matched->p4FromVtx(vtxPos,matched->finalEnergy).Pt();
      outRecoEtaSC = matched->SC.eta;
      outR9 = matched->SC.r9;

      if(!etaSelectPhoton(matched)){
	outputTree->Fill(); continue;
      }
      passEta=true;
      
      if(!PhotonID->getPreSelection(matched,nVtx,rho,0)){
	outputTree->Fill(); continue;
      }
      passPreselection=true;

      if(!PhotonID->getIdCiCPF(matched,nVtx,rho,0)){
	outputTree->Fill(); continue;
      }
      passCiC=true;
      
      outputTree->Fill();
    }

}

void HggEfficiencyMap::runMMG(){
  TVector3 vtxPos(vtxX[0],vtxY[0],vtxZ[0]); 
  
  //loop over muons
  for(int iMu=0;iMu<nMu_;iMu++){
    VecbosMu& mu1 = Muons_->at(iMu);
    if(mu1.pt<10) continue;
    if(!HggMuonID::passMuMuGammaID(mu1)) continue;

    for(int jMu=iMu+1;jMu<nMu_;jMu++){
      VecbosMu& mu2 = Muons_->at(jMu);
      if(mu2.pt<10) continue;
      if(!HggMuonID::passMuMuGammaID(mu2)) continue;
      if(mu1.pt<20 && mu2.pt<20) continue; //require 20/10

      for(int iPho=0;iPho<nPho_;iPho++){
	VecbosPho& pho = Photons_->at(iPho);

	TLorentzVector p4Mu1; p4Mu1.SetPtEtaPhiM(mu1.pt,mu1.eta,mu1.phi,0.106);
	TLorentzVector p4Mu2; p4Mu2.SetPtEtaPhiM(mu2.pt,mu2.eta,mu2.phi,0.106);
	TLorentzVector p4Pho = pho.p4FromVtx(vtxPos,pho.finalEnergy);
	float Mmumu = (p4Mu1+p4Mu2).M();

	float M = (p4Mu1+p4Mu2+p4Pho).M();
	if(Mmumu>75 || M<88 || M>94) continue;
	this->clear();
	outPt = -1;
	outEta = 0;
	outEtaSC = 0;
	
	outRecoPt=p4Pho.Pt();
	outRecoEta=p4Pho.Eta();;
	outRecoEtaSC=pho.SC.eta;
	outPhi = pho.SC.phi;
	outR9=pho.SC.r9;
	if(!etaSelectPhoton(&pho)){
	  outputTree->Fill(); continue;
	}
	passEta=true;
	
	if(!PhotonID->getPreSelection(&pho,nVtx,rho,0)){
	  outputTree->Fill(); continue;
	}
	passPreselection=true;
	
	if(!PhotonID->getIdCiCPF(&pho,nVtx,rho,0)){
	  outputTree->Fill(); continue;
	}
	passCiC=true;
	
	outputTree->Fill();	
      }
    }
  }
}

void HggEfficiencyMap::runZee(){
  TVector3 vtxPos(vtxX[0],vtxY[0],vtxZ[0]); 
  for(int i=0;i<nPho_;i++){
    VecbosPho& pho1 = Photons_->at(i);
    //tag electron
    if(!PhotonID->getIdCiCPF(&pho1,nVtx,rho,0)) continue;
    TLorentzVector p4ele1 = pho1.p4FromVtx(vtxPos,pho1.finalEnergy);
    if(p4ele1.Pt() < 30) continue;
    for(int j=0;j<nPho_;j++){
      if(j==i) continue; //find another photon consistent with the Z peak
      VecbosPho& pho2 = Photons_->at(j);
      TLorentzVector p4ele2 = pho2.p4FromVtx(vtxPos,pho2.finalEnergy);
      float M = (p4ele1+p4ele2).M();

      if(M<88 || M > 94) continue;

      //OK, photon j is consistent with the Z peak
      this->clear();

	outPt = -1;
	outEta = 0;
	outEtaSC = 0;
	
	outRecoPt=p4ele2.Pt();
	outRecoEta=p4ele2.Eta();;
	outRecoEtaSC=pho2.SC.eta;
	outPhi = pho2.SC.phi;
	outR9=pho2.SC.r9;
	if(!etaSelectPhoton(&pho2)){
	  outputTree->Fill(); continue;
	}
	passEta=true;
	
	if(!PhotonID->getPreSelection(&pho2,nVtx,rho,0)){
	  outputTree->Fill(); continue;
	}
	passPreselection=true;
	
	if(!PhotonID->getIdCiCPF(&pho2,nVtx,rho,0)){
	  outputTree->Fill(); continue;
	}
	passCiC=true;
	
	outputTree->Fill();	
    }
  }  
}

bool HggEfficiencyMap::etaSelectPhoton(VecbosPho* pho){
  //apply kinematic photon selection
  if(fabs(pho->SC.eta) > 2.5) return false;
  if(fabs(pho->SC.eta) > 1.4442 && fabs(pho->SC.eta) < 1.566) return false;
  
  return true;
}


void HggEfficiencyMap::smearPhoton(VecbosPho* pho,int smearShift){
  TRandom3 rng(0);
  float smear = pho->dEoE + smearShift*pho->dEoEErr;
  if(smear < 0) smear = 0;
  if(smearShift<-100) smear = 0;
  float rand=0;
  if(smear > 0) rand = rng.Gaus(0,smear);
  if(rand < -1) rand=-1;
  if(rand > 1E3) rand = 1E3;
  pho->finalEnergy = pho->scaledEnergy*(1+rand);
  pho->finalEnergyError = pho->scaledEnergyError*(1+rand);
}

void HggEfficiencyMap::clear(){
    passMatch=false;
    passEta=false;
    passPreselection=false;
    passCiC=false;

}

void HggEfficiencyMap::setBranchAddresses(){
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

 fChain->SetBranchAddress("nMu",&nMu_);
 fChain->SetBranchAddress("Muons",&Muons_);

 fChain->SetBranchAddress("nJet",&nJet_);
 fChain->SetBranchAddress("Jets",&Jets_);

 fChain->SetBranchAddress("nGenHiggs",&nGenHiggs);
 fChain->SetBranchAddress("GenHiggs",&GenHiggs);
 
 fChain->SetBranchAddress("nGenPho",&nGenPho);
 fChain->SetBranchAddress("GenPhotons",&GenPhotons);

 fChain->SetBranchAddress("nGenEle",&nGenEle);
 fChain->SetBranchAddress("GenElectrons",&GenElectrons);

fChain->SetBranchAddress("nPU",&inPU);

 fChain->SetBranchAddress("PFMET",&pfMet);
 fChain->SetBranchAddress("PFMETPhi",&pfMetPhi);

}


VecbosPho* HggEfficiencyMap::findGenMatchPhoton(VecbosGen* gen){
  int index = gen->index;

  for(int iPho=0;iPho<nPho_;iPho++){
    VecbosPho* pho = &(Photons_->at(iPho));
    if(pho->genMatch.index == index) return pho;
  }
  return 0;
}

VecbosPho* HggEfficiencyMap::findGenMatchElectron(VecbosGen* gen){
  //need to match a RECO photon to the gen level electron
  const float maxDR = 0.2;
  float dEoEBest = 9999;
  int indexPho = -1;

  for(int iPho=0;iPho<nPho_;iPho++){
    VecbosPho* pho = &(Photons_->at(iPho));
    if(DeltaR(pho->SC.eta,gen->eta,pho->SC.phi,gen->phi) > maxDR) continue;
    
    float dEoE = fabs(pho->finalEnergy-gen->energy)/gen->energy;
    if(dEoE > 1.) continue;
    if(dEoE < dEoEBest){
      dEoEBest = dEoE;
      indexPho = iPho;
    }
  }  
  if(indexPho==-1) return 0;
  return &(Photons_->at(indexPho));
}
