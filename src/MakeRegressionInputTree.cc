#include "MakeRegressionInputTree.hh"



MakeRegressionInputTree::MakeRegressionInputTree(std::vector<std::string> fNames, std::string treeName,
			std::string outputFile):
  BaseSelector(fNames,treeName,outputFile) 
{
  setDoFill(false); //we will manage the filling of the output tree
  }

int MakeRegressionInputTree::init() {
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("nPho",1);
  fChain->SetBranchStatus("nEle",1);
  fChain->SetBranchStatus("Photons.*",1);
  fChain->SetBranchStatus("Electrons.*",1);
  fChain->SetBranchStatus("Photons.correctedEnergy",1);
  fChain->SetBranchStatus("Photons.correctedEnergyError",1);
  fChain->SetBranchStatus("Photons.SC.eta",1);
  fChain->SetBranchStatus("Photons.SC.r9",1);
  fChain->SetBranchStatus("Photons.phi",1);
  fChain->SetBranchStatus("Photons.eta",1);
  fChain->SetBranchStatus("Photons.genMatch.index",1);
  fChain->SetBranchStatus("Photons.dr03HcalTowerSumEtCone",1);
  fChain->SetBranchStatus("Photons.dr03TrackIso[100]",1);
  fChain->SetBranchStatus("Photons.dr02ChargedHadronPFIso[100]",1);
  fChain->SetBranchStatus("Photons.HoverE",1);
  fChain->SetBranchStatus("Photons.SC.sigmaIEtaIEta",1);
  
  fChain->SetBranchStatus("Photons.dr03ChargedHadronPFIso[100]",1);
  fChain->SetBranchStatus("Photons.dr03PhotonPFIso",1);
  
  fChain->SetBranchStatus("nPU",1);
  fChain->SetBranchStatus("nVtx",1);
  fChain->SetBranchStatus("rho",1);
  
  fChain->SetBranchStatus("photonMatchedElectron",1);

  return 0;
}

void MakeRegressionInputTree::processEntry(Long64_t iEntry) {
  pu=inPU;

  float pass_pt[2]={-1,-1}; //store the pts of the passing photons
  int   pass_index[2]={-1,-1};
  for(int i=0;i<nPho_;i++) {
    auto photon = &(Photons_->at(i));
    pt = photon->correctedEnergy/cosh(photon->eta);
    passWP90_id  = photonID.passID (*photon,StandardPhotonID::kLoose);
    passWP90_iso = photonID.passIso(*photon,StandardPhotonID::kLoose);
    passWP90 = passWP90_id && passWP90_iso;

    se = photon->correctedEnergyError/photon->correctedEnergy;
    etaSC = photon->SC.eta;
    r9 = photon->SC.r9;
    phi = photon->phi;
    pt  = photon->correctedEnergy/cosh(photon->eta);

    realPho = (photon->genMatch.index != -1);
    //check if its a real electron
    realEle = false;
    for(int iEle=0;iEle<nEle_;iEle++){
      if( Electrons_->at(iEle).SC.index = photon->SC.index){
	if( Electrons_->at(iEle).genMatch.index ) realEle=true;
      }
    }
    
    pfPhoton = photon->dr03PhotonPFIso;
    pfChargedGood = photon->dr03ChargedHadronPFIso[0];

    HE = photon->HoverE;
    sieie = photon->SC.sigmaIEtaIEta;
    sieip = photon->SC.sigmaIEtaIPhi;
    sipip = photon->SC.sigmaIPhiIPhi;

    e3x3 = photon->SC.e3x3;
    e5x5 = photon->SC.e5x5;
    rawE = photon->SC.rawE;
    etaWidth = photon->SC.etaWidth;
    phiWidth = photon->SC.phiWidth;
    
    nBC = photon->SC.nBCs;
    energyBC = photon->SC.BCSeed.energy;
    etaBC = photon->SC.BCSeed.eta;
    phiBC = photon->SC.BCSeed.phi;
    
    eMax = photon->SC.eMax;
    e2nd = photon->SC.e2nd;
    eTop = photon->SC.eTop;
    eBottom = photon->SC.eBottom;
    eLeft = photon->SC.eLeft;
    eRight = photon->SC.eRight;
    
    e2x5Max = photon->SC.e2x5Max;
    e2x5Top = photon->SC.e2x5Top;
    e2x5Bottom = photon->SC.e2x5Bottom;
    e2x5Left = photon->SC.e2x5Left;
    e2x5Right = photon->SC.e2x5Right;
    
    electronMatch = photonMatchedElectron[i];

    outTree->Fill();
  }
}

void MakeRegressionInputTree::clear() {}

void MakeRegressionInputTree::setupOutputTree() {
  outTree = new TTree("RegressionInputTree","");

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
  outTree->Branch("passWP90",&passWP90);
  outTree->Branch("passWP90_id",&passWP90_id);
  outTree->Branch("passWP90_iso",&passWP90_iso);
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
  outTree->Branch("nPU",&pu,"nPU/I");
}
