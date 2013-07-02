// std includes

#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TH1D.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "ReadConfig.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "VecbosEGObject.hh"
#include "HggReducer.hh"
#include "../src/HggPhysUtils.cc"
#include "assert.h"

#define debugReducer 0

HggReducer::HggReducer(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
  vertexCFG = "";
  energyScaleCFG = "";
  energySmearCFG = "";
  mcScalingCFG = "";
  minPhoSel = 2;
}

HggReducer::HggReducer(TTree *tree, string json, bool goodRunLS, bool isData,int mod) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;
  vertexCFG = "";
  energyScaleCFG = "";
  energySmearCFG = "";
  mcScalingCFG = "";
  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
  minPhoSel = 2;
}

void HggReducer::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void HggReducer::SetWeight(double weight) {
  _weight = weight;
}

struct less_than_pt_photon{
  inline bool operator() (const VecbosPho& p1, const VecbosPho& p2){ 
    return p1.finalEnergy/cosh(p1.eta) < p2.finalEnergy/cosh(p2.eta);
  }
};

void HggReducer::Loop(string outFileName, int start, int stop) {
  if(fChain == 0){
    cout << "fChain Not defined! QUITTING" << endl;
    return;
  }
  this->init(outFileName);    


  //do initializations (setup output tree, etc.)

  if(debugReducer) cout << "Doing Initialization ... " << flush;
  if(debugReducer) cout << "Done" <<endl;
  
  TRandom3 rng(0);

  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // setup for triggers:
  if(debugReducer) cout << "setting Triggers: " << endl;
  const int nTrigs = triggerNames.size();
  vector<vector<string> > masks;
  for(int iTrig=0;iTrig<nTrigs;iTrig++){
    if(debugReducer) cout << ">> " << triggerNames.at(iTrig) << endl; 
    std::vector<string> tmp;
    tmp.push_back(triggerNames.at(iTrig));// setup the trigger masks for Vecbos
    masks.push_back(tmp);
  }

  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Starting with Entry: " << start << endl;
  cout << "Number of entries = " << stop << endl;
  
  Long64_t jentry = start-1;
  while(fChain->GetEntry(++jentry)){
    if(jentry == stop) break;
    if (jentry%500 == 0) cout << ">>> Processing event # " << jentry << endl;

    runNumberO=runNumber;
    evtNumberO=eventNumber;
    lumiBlockO=lumiBlock;

    //filters:
    if(nTrack > 5000 || nPV > 160) {
      if(debugReducer) cout << "dropping event: too many tracks/PVs: " << nTrack << "  " << nPV << endl; 
      continue;
    }

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      // hadronic reload trigger masks and get HLT decisions
      for(int iTrig=0;iTrig<nTrigs;iTrig++){
	setRequiredTriggers(masks.at(iTrig));
	reloadTriggerMask(true);
	triggerBits[iTrig] = hasPassedHLT();
      }
    }

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      if(debugReducer) cout << "Dropping: " << runNumber << ":" << lumiBlock << endl;
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }    

    eeBadScFilterFlag         = METFlags & (1 << 8);
    hcalLaserEventFilterFlag  = METFlags & (1 << 7);
    HBHENoiseFilterResultFlag = METFlags & (1 << 6);
    isNotDeadEcalCluster      = METFlags & (1 << 5);
    trackerFailureFilterFlag  = METFlags & (1 << 4);
    CSCHaloFilterFlag         = METFlags & (1 << 3);
    drDead                    = METFlags & (1 << 2); 
    drBoundary                = METFlags & (1 << 1);
    ECALTPFilterFlag          = METFlags & (1 << 0);
    
    this->clearAll();
    this->fillVertexInfo();

    this->fillMuons();
    this->fillElectrons();
    this->fillJets();
    // setup the photons
    std::vector<VecbosPho> tmpPhotons; //temporarily hold a collection of VecbosPhos
    int nSelectedPho=0; // number of selected photons
    if(debugReducer) cout << ">> " << nPho << " Photons" << endl;

    for(int iPho=0;iPho<nPho;iPho++){
      VecbosPho pho(this,iPho);
      //if(!_isData) scaler->ScalePhoton(pho);
      if(debugReducer) cout << iPho << ": energy=" <<  pho.energy << "  eta=" << pho.SC.eta << endl;
      if(pho.energy/cosh(pho.eta) < 10) continue;
      if(fabs(pho.SC.eta) > 3.) continue;
      //if(!this->passPreselection(&pho)) continue;
      if(debugReducer) cout << "pass" << endl;

      //DO ENERGY CORRECTION
      corrector->getPhotonEnergyCorrection(pho,false);
      if(debugReducer) cout << "Corrected Photon: E=" << pho.energy << "  CorE=" << pho.correctedEnergy << "  eta=" 
			    << pho.SC.eta << "  r9=" << pho.SC.r9 << endl;
      


      if(_isData){ // get the energy scale for data
	if(debugReducer) cout << "Doing Energy Scale ... " << flush;
	std::pair<float,float> dE = energyScale->getDEoE(pho,runNumber);
	if(debugReducer) cout << "Done" << endl;;	
	pho.dEoE    = dE.first;
	pho.dEoEErr = 0; 
	pho.scaledEnergy = pho.correctedEnergy*(pho.dEoE);
	pho.scaledEnergyError = pho.correctedEnergyError*((pho.dEoE+pho.dEoEErr));
	if(debugReducer) cout << pho.dEoE <<"   " << pho.scaledEnergy << endl;
      
      }


      if(!_isData){ // monte carlo, get energy smearing and scale error
	//first get the scale error
	pho.scaledEnergy = pho.correctedEnergy;
	pho.scaledEnergyError = pho.scaledEnergy*energySmear->getMCScaleErr(pho,applyScaleSmear);
	
	//now get the energy smearing
	std::pair<float,float> dE = energySmear->getDEoE(pho,applyScaleSmear);	
	pho.dEoE    = dE.first;
	pho.dEoEErr = dE.second;  // from these, we can generate the energy smearing later (its just a gaussian)
      }

      pho.finalEnergy = pho.scaledEnergy;
      pho.finalEnergyError = pho.scaledEnergyError;

      nSelectedPho++;
      
      // ADD things to the collections
      Photons_.push_back(pho);
      nPho_++;
      if(debugReducer) cout << "aa " << nPho_ << endl;
      //Photons_.push_back(pho);
      //nPho_++;
      /*
      VecbosSC sc = pho.SC;
      SuperClusters_.push_back(sc);
      nSC_++;
      VecbosPFSC pfsc = pho.PFSC;
      PFSuperClusters_.push_back(pfsc);      
      nPFSC++;
      //for(int i=0;i<4;i++){ // push back up to 4 BCs for each SC
      //	if(scSt.BCs[i] > -1) BasicClusters_.push_back(pho.SC.basicClusters[i].getStruct());
      //}
      if(pho.hasConv()){ //there is a valid conversion
	Conversions_.push_back(pho.conversions.at(0));
	nConv_++;
      }
      */
      //end of photon loop
    }//for(int iPho=0; ...

    //now sort the photons by pT
    std::sort(Photons_.begin(),Photons_.end(),less_than_pt_photon());


    
    if(nSelectedPho<minPhoSel){
      if(debugReducer) cout << "Fewer than " << minPhoSel << " selected photons" << endl;
      continue; // skip the event if there are fewer than 2 photons
    }
    
    if(debugReducer) cout << "More than 1 selected photons" << endl;
    //do the vertexing TMVA
    std::vector<VecbosPho>::iterator iPho1;
    std::vector<VecbosPho>::iterator iPho2;
    if(debugReducer) cout << "aa" << endl;
    nPair_=0;
    for(iPho1 = Photons_.begin(); iPho1 != Photons_.end(); iPho1++){
      if(debugReducer) cout << nPair_ << endl;
      //VecbosPho pho1 = *iPho1;
      for(iPho2 = iPho1+1; iPho2 != Photons_.end(); iPho2++){
	if(debugReducer) cout << ">> " << nPair_ << endl;
	if(iPho2==iPho1) continue;
	if(debugReducer) cout << ">>>>" << endl;
	//VecbosPho pho2 = *iPho2;
	if(debugReducer) cout << "Doing Vertexing ... " << endl;
	float perEvt=-1;
	vector<pair<int,float> > vtxPair = vertexer->vertex_tmva(&*iPho1,&*iPho2,perEvt);

	const int nTop=3;
	pair<int,float> top[nTop];
	for(int i=0;i<nTop;i++) top[i] = pair<int,float>(0,-1);
	vector<pair<int,float> >::const_iterator vtxIt;
	for(vtxIt = vtxPair.begin(); vtxIt != vtxPair.end(); vtxIt++){
	  pair<int,float> tmp = *vtxIt;
	  for(int i=0;i<nTop;i++){
	    if(tmp.second > top[i].second){
	      swap(tmp,top[i]);
	    }
	  }
	}
	ggVerticesPhotonIndices.push_back(pair<int,int>(iPho1->index,iPho2->index) );
	ggVerticesVertexIndex01.push_back(top[0]);
	ggVerticesVertexIndex02.push_back(top[1]);
	ggVerticesVertexIndex03.push_back(top[2]);
	assert(top[0].second >= top[1].second); //make sure the vertices are in descending order

	ggVerticesPerEvtMVA.push_back(perEvt);
	nPair_++;
      }//for(iPho2 = iPho1+1; iPho2 != Photons_.end(); iPho2++)
    }//for(iPho1 = Photons_.begin(); iPho1 != Photons_.end(); iPho1++)
    
    // fill the isolation variables for the photons
    if(debugReducer) cout << "Filling IsoVars" << endl;
    for(int iPho = 0; iPho < nPho_; iPho++){
      Photons_.at(iPho).nPV=nVtx;
      for(int i=0;i<nVtx;i++){
	float thisIsoDR02 = computeTrackIso(iPho,i,1e-5,0.2,0.02,0.0,1.0,0.1);
	float thisIsoDR03 = computeTrackIso(iPho,i,1e-5,0.3,0.2,0.0,1.0,0.1);
	float thisIsoDR04 = computeTrackIso(iPho,i,1e-5,0.4,0.2,0.0,1.0,0.1);
	Photons_.at(iPho).dr02TrackIso[i]=thisIsoDR02;
	Photons_.at(iPho).dr03TrackIso[i]=thisIsoDR03;
	Photons_.at(iPho).dr04TrackIso[i]=thisIsoDR04;
      }
    } // end of photon isolation loop

    //fill remaining collections
    if(!_isData) this->fillGeneratorInfo();
    this->matchPhotonsElectrons();
   //cout << "Filling Collection" << endl;
    outTree->Fill();
  } // end of main loop

  cout << "Writing Tree:" << endl;

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();

  //write output
 file->Close();
}


float HggReducer::computeTrackIso(int iPho, int iVtx,
				  float ptMin,
				  float outerCone,
				  float innerCone,
				  float etaStripHalfWidth,
				  float dzMax,
				  float dxyMax){

  float vX=0,vY=0,vZ=0;
  if(iVtx >= 0 && iVtx < nPV){
    vX = vtxX[iVtx];
    vY = vtxY[iVtx];
    vZ = vtxZ[iVtx];
  }

  double SumTrackPt=0;
  for(int iTrack=0;iTrack<nTrack;iTrack++){
    double ptTrack = sqrt( pxTrack[iTrack]*pxTrack[iTrack] + pyTrack[iTrack]*pyTrack[iTrack] );
    if(ptTrack < ptMin) continue;
    double dZ = fabs( (trackVzTrack[iTrack] - vZ) - 
		      ( (trackVxTrack[iTrack] - vX)*pxTrack[iTrack] 
			+ (trackVyTrack[iTrack] - vY)*pyTrack[iTrack])/ptTrack*pzTrack[iTrack]/ptTrack);
    if(dZ > dzMax) continue;
    double dXY = ( (trackVyTrack[iTrack] - vY)*pyTrack[iTrack] - (trackVxTrack[iTrack] - vX)*pxTrack[iTrack] )/ptTrack;
    if( fabs(dXY) > dxyMax ) continue;
    TVector3 trackP(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
    if(trackP.Pt()<1e-4) continue;
    assert(trackP.Pt()>0);
    assert(trackP.Phi() > -4 && trackP.Phi()<4);
    assert(phiPho[iPho] > -4 && phiPho[iPho]<4);
    double dEta = fabs(etaPho[iPho] - trackP.Eta());
    double dR   = DeltaR<float>(etaPho[iPho],trackP.Eta(),phiPho[iPho],trackP.Phi());
    if( dR < outerCone && dR >= innerCone && dEta >= etaStripHalfWidth) SumTrackPt+=ptTrack;
  }
  return SumTrackPt;
}

void HggReducer::matchPhotonsElectrons(){
  vector<VecbosPho>::const_iterator phoIt;
  for(phoIt = Photons_.begin(); phoIt != Photons_.end();phoIt++){
    bool match = false;
    for(int i=0;i<nEle;i++){ 
      //copied from  ConversionTools::hasMatchedPromptElectron
      if(superClusterIndexEle[i] != phoIt->SC.index) continue;
      if(hasMatchedConversionEle[i]) continue;
      if(gsfTrackIndexEle[i]<0 || gsfTrackIndexEle[i] >=nGsfTrack) continue;
      if(expInnerLayersGsfTrack[gsfTrackIndexEle[i]]>0) continue;
      match = true;
      break;	
    }
    photonMatchedElectron[phoIt-Photons_.begin()] = match;
  }
}

void HggReducer::clearAll(){
//clear collections

  nPho_=0;
  Photons_.clear();
  
  nMu_=0;
  Muons_.clear();
  nEle_=0;
  Electrons_.clear();
  nJet_=0;
  Jets_.clear();
  pileupBunchX->clear();
  pileupNInteraction->clear();
  
  ggVerticesPhotonIndices.clear();
  ggVerticesVertexIndex01.clear();
  ggVerticesVertexIndex02.clear();
  ggVerticesVertexIndex03.clear();
  ggVerticesPerEvtMVA.clear();

  nGenPho=0;
  nGenMu =0;
  nGenEle=0;
  nGenHiggs=0;

  GenPhotons.clear();
  GenMuons.clear();
  GenElectrons.clear();
  GenHiggs.clear();

}

void HggReducer::init(string outputFileName){
  //define the tree
  outTree = new TTree("HggReduce","HggReduce");

 //setup the trigger objects
 triggerBits = new int[triggerNames.size()];

 pileupBunchX = new std::vector<short>; pileupBunchX->clear();
 pileupNInteraction = new std::vector<short>; pileupNInteraction->clear();

 this->setOutputBranches(); // initialize the branches of the output tree

 ReadConfig cfg;
 try{
   cfg.read(config);
 }catch(exception &e){
   std::cout << "HggReducer" << std::endl;
   std::cout << "config file: " << config <<std::endl;
   throw e;
   
 }
 string EnergyScaleCFG  = cfg.getParameter("EnergyScaleCFG");
 string EnergySmearCFG  = cfg.getParameter("EnergySmearCFG");
 string MCScalingCFG    = cfg.getParameter("ScalingFile");
 string sScaleSmear     = cfg.getParameter("ScaleSmear");
 string sMinPhoSel      = cfg.getParameter("MinPreselPhotons");
 string puWeight        = cfg.getParameter("pileupReweight");
 triggerNames  = cfg.getTokens("Triggers",",");

 applyScaleSmear = atoi(sScaleSmear.c_str());
 if(sMinPhoSel.compare("")!=0) minPhoSel = atoi(sMinPhoSel.c_str());

 pileupWeight=1;
 if(puWeight.compare("")!=0){
   cout << "Opening Pileup Weight File: " << puWeight << endl;
   pileupWeightFile = new TFile(puWeight.c_str());
   pileupWeightHist = (TH1F*)pileupWeightFile->Get("pileupReWeight");
 }

 correctJets = cfg.getParameter("doJetCorrection").compare("yes")==0;

 cout << "Config Parameters:" << endl
      << "EnergyScale: " << EnergyScaleCFG << endl
      << "EnergySmear: " << EnergySmearCFG << endl
      << "MC Scaling:  " << MCScalingCFG << endl
      << "ApplyScaleSmear: "  << sScaleSmear << endl
      << "Requiring " << minPhoSel << " Photons" <<endl
      << "Doing Jet Corrections: " << correctJets << endl;

 
 vertexer  = new HggVertexing(this);
 vertexer->setConfigFile(config);
 vertexer->useConversions();
 vertexer->saveInputs(outTree);
 vertexer->init();
 
 corrector = new HggEGEnergyCorrector(this,config,_isData);
 elecorrector = new HggEGEnergyCorrector(this,config,_isData);
 elecorrector->useElectronWeights();
 energyScale = new HggEnergyScale(EnergyScaleCFG);
 energySmear = new HggEnergyScale(EnergySmearCFG);
 
 if(correctJets) jetCorr = new VecbosJetCorrector(cfg);
 
 
 if(!_isData) scaler = new HggScaling(MCScalingCFG);
}


void HggReducer::fillMuons(){
  nMu_=0;
  for(int iMuon = 0; iMuon<nMuon;iMuon++){
    VecbosMu mu(this,iMuon);
    if(mu.index==-1) continue;
    Muons_.push_back(mu);
    nMu_++;
  }
}
void HggReducer::fillElectrons(){
  nEle_=0;
  for(int iEle = 0; iEle<nEle;iEle++){
    VecbosEle ele(this,iEle);
    if(ele.index==-1) continue;
    elecorrector->getElectronEnergyCorrection(ele);
    Electrons_.push_back(ele);
    nEle_++;
  }
}

void HggReducer::fillJets(){
  const float minPt = 5;
  nJet_=0;
  for(int iJet = 0; iJet<nAK5PFPUcorrJet;iJet++){
    VecbosJet jet(this,iJet,VecbosJet::PFPUcorr);
    if(jet.index==-1) continue;
    if(jet.pt < minPt) continue;
    if(correctJets) jetCorr->CorrectJet(jet,rhoJetsFastJet);
    Jets_.push_back(jet);
    nJet_++;
  }
  caloMet =  TMath::Sqrt(TMath::Power(pxMet[0],2)+TMath::Power(pyMet[0],2));
  caloMetPhi =  phiMet[0];
  pfMet = TMath::Sqrt(TMath::Power(pxPFMet[0],2)+TMath::Power(pyPFMet[0],2));
  pfMetPhi = phiPFMet[0];
  tcMet =  TMath::Sqrt(TMath::Power(pxTCMet[0],2)+TMath::Power(pyTCMet[0],2));
  tcMetPhi = phiTCMet[0];

  return;
  const int nJetCat=4;
  const float maxJetEta[nJetCat] = {2.5,2.75,3,4.7};
  const float betaStarSlope[nJetCat] = {0.2,0.3,999,999};
  const float rmsCut[nJetCat] = {0.06,0.05,0.05,0.055};

  //nJets=0;
  for(int iJ=0;iJ<nAK5PFPUcorrJet; iJ++){
    int jetCat=0;
    for(; jetCat<4; jetCat++){ // get the jet category
      if(fabs(etaAK5PFPUcorrJet[iJ]) <maxJetEta[jetCat]) break;
    }
    if(jetCat>3) continue; //if its >4.7, reject the jet
    float pT = TMath::Sqrt(TMath::Power(pxAK5PFPUcorrJet[iJ],2)+TMath::Power(pyAK5PFPUcorrJet[iJ],2));
    if(pT < minPt) continue;
    if(betastarIdMvaAK5PFPUcorrJet[iJ] > betaStarSlope[jetCat]*TMath::Log(nPV)-0.64) continue;
    if(rmsCandsHandAK5PFPUcorrJet[iJ] > rmsCut[jetCat]) continue; //jet ID variables

    /*    
    ptJet[nJets] = pT;
    etaJet[nJets] = etaAK5PFPUcorrJet[iJ];
    phiJet[nJets] = phiAK5PFPUcorrJet[iJ];
    energyJet[nJets] = energyAK5PFPUcorrJet[iJ];
    nJets++;
    */
  }
}
void HggReducer::fillVertexInfo(){
  nVtx = nPV;
  for(int i=0;i<nPV;i++){
    vtxX[i] = PVxPV[i];             
    vtxY[i] = PVyPV[i];             
    vtxZ[i] = PVzPV[i];             
    vtxChi2[i] = chi2PV[i];          
    vtxNdof[i] = ndofPV[i];          
    vtxNormalizedChi2[i] = normChi2PV[i];
    vtxTrackSize[i] = trackSizePV[i];       
    vtxIsFake[i] = isFakePV[i];          
    vtxIsValid[i] = isValidPV[i];         
  }
  rho=rhoJetsFastJet;  // rho from kt6PFJets
}

void HggReducer::fillGeneratorInfo(){
  //GENERATOR information
  nGenPho=0;
  nGenMu =0;
  nGenEle=0;
  nGenOthers=0;
  nGenHiggs=0;

  for(int iGen=0;iGen<nMc;iGen++){
    VecbosGen part(this,iGen);

    switch(abs(idMc[iGen])){
    case 25:  //higgs
      GenHiggs.push_back(part); break;
    case 22: //Photons
      GenPhotons.push_back(part); break;
    case 13: //Muons
      GenMuons.push_back(part); break;
    case 11: //Electrons
      GenElectrons.push_back(part); break;
    Default:
      GenOthers.push_back(part); break;
    }
    nGenHiggs = GenHiggs.size();
    nGenPho   = GenPhotons.size();
    nGenEle   = GenElectrons.size();
    nGenMu    = GenMuons.size();
    nGenOthers= GenOthers.size();
  }    
  procID = 0; //genProcessId;
  qScale;
  nPu = nPU[1];
  pileupWeight = pileupWeightHist->GetBinContent( pileupWeightHist->FindFixBin(nPu) );

  //match the reco objects to gen objects
  
  PhoCollection::iterator phoIt;
  for(phoIt = Photons_.begin(); phoIt !=Photons_.end(); phoIt++)
    phoIt->doGenMatch(this);
  MuCollection::iterator muIt;
  for(muIt = Muons_.begin(); muIt !=Muons_.end(); muIt++)
    muIt->doGenMatch(this);
  EleCollection::iterator eleIt;
  for(eleIt = Electrons_.begin(); eleIt !=Electrons_.end(); eleIt++)
    eleIt->doGenMatch(this);

}

void HggReducer::setOutputBranches(){

  //Event info
  outTree->Branch("lumiBlock",&lumiBlockO,"lumiBlock/I");
 outTree->Branch("runNumber",&runNumberO,"runNumber/I");
outTree->Branch("evtNumber",&evtNumberO,"evtNumber/I");
 outTree->Branch("bunchX",&bunchX,"bunchX/I");
 outTree->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
 outTree->Branch("evtTime",&evtTime,"evtTime/I");
 outTree->Branch("isRealData",&_isData,"isRealData/I");

 // MET Flags
 outTree->Branch("eeBadScFilterFlag",&eeBadScFilterFlag,"eeBadScFilterFlag/B");
 outTree->Branch("hcalLaserEventFilterFlag",&hcalLaserEventFilterFlag,"hcalLaserEventFilterFlag/B");
 outTree->Branch("HBHENoiseFilterResultFlag",&HBHENoiseFilterResultFlag,"HBHENoiseFilterResultFlag/B");
 outTree->Branch("isNotDeadEcalCluster",&isNotDeadEcalCluster,"isNotDeadEcalCluster/B");
 outTree->Branch("trackerFailureFilterFlag",&trackerFailureFilterFlag,"trackerFailureFilterFlag/B");
 outTree->Branch("CSCHaloFilterFlag",&CSCHaloFilterFlag,"CSCHaloFilterFlag/B");
 outTree->Branch("drDead",&drDead,"drDead/B");
 outTree->Branch("drBoundary",&drBoundary,"drBoundary/B");
 outTree->Branch("ECALTPFilterFlag",&ECALTPFilterFlag,"ECALTPFilterFlag/B");

 ///information for the vertex
 outTree->Branch("nVtx",&nVtx,"nVtx/I");
 outTree->Branch("vtxX",vtxX,"vtxX[nVtx]/F");
 outTree->Branch("vtxY",vtxY,"vtxY[nVtx]/F");
 outTree->Branch("vtxZ",vtxZ,"vtxZ[nVtx]/F");
 outTree->Branch("vtxChi2",vtxChi2,"vtxChi2[nVtx]/F");
 outTree->Branch("vtxNdof",vtxNdof,"vtxNdof[nVtx]/F");
 outTree->Branch("vtxNormalizedChi2",vtxNormalizedChi2,"vtxNormalizedChi2[nVtx]/F");
 outTree->Branch("vtxTrackSize",vtxTrackSize,"vtxTrackSize[nVtx]/I");
 outTree->Branch("vtxIsFake",vtxIsFake,"vtxIsFake[nVtx]/I");
 outTree->Branch("vtxIsValid",vtxIsValid,"vtxIsValid[nVtx]/I");
 
 //physics declared -- should be set by the JSON, but can't hurt
 outTree->Branch("phyDeclared",&phyDeclared,"phyDeclared/I");
 
 
 outTree->Branch("rho", &rho,"rho/F");
 outTree->Branch("rhoEtaMax44", &rhoEtaMax44,"rhoEtaMax44/F");
 
 
 
 outTree->Branch("pileupBunchX","std::vector<short>", &pileupBunchX);
 outTree->Branch("pileupNInteraction","std::vector<short>", &pileupNInteraction);
 outTree->Branch("pileupTrueNumInterations",&pileupTrueNumInterations);
 

  //trigger -- here we depart a bit from Yong's original code
 for(int i=0;i<triggerNames.size();i++){
   outTree->Branch(triggerNames.at(i).c_str(), &(triggerBits[i]), Form("%s/I",triggerNames.at(i).c_str()) );  // this will produce 1 int per trigger in the output tree
 }
 

 //objects
 /*
 outTree->Branch("nSC",&nSC_);
 outTree->Branch("nPFSC",&nPFSC_);
 outTree->Branch("nConv",&nConv_);
 outTree->Branch("SuperClusters",&SuperClusters_);                          
 outTree->Branch("PFSuperClusters",&PFSuperClusters_);                          
 //outTree->Branch("BasicClusters",&BasicClusters_);                          
 outTree->Branch("Conversions",&Conversions_);
 */
 outTree->Branch("nPho",&nPho_);
 outTree->Branch("nPair",&nPair_);
 outTree->Branch("Photons",&Photons_);
 outTree->Branch("ggVerticesPhotonIndices",&ggVerticesPhotonIndices);
 outTree->Branch("ggVerticesVertexIndex",&ggVerticesVertexIndex01);
 outTree->Branch("ggVerticesVertexIndex2nd",&ggVerticesVertexIndex02);
 outTree->Branch("ggVerticesVertexIndex3rd",&ggVerticesVertexIndex03);
 outTree->Branch("ggVerticesPerEvtMVA",&ggVerticesPerEvtMVA);
 outTree->Branch("photonMatchedElectron",photonMatchedElectron,"photonMatchedElectron[nPho]/O");

 outTree->Branch("nMu",&nMu_,"nMu/I");
 outTree->Branch("Muons",&Muons_);

 outTree->Branch("nEle",&nEle_,"nEle/I");
 outTree->Branch("Electrons",&Electrons_);

 outTree->Branch("nJet",&nJet_,"nJet/I");
 outTree->Branch("Jets",&Jets_);

 outTree->Branch("CaloMET",&caloMet,"CaloMET");
 outTree->Branch("CaloMETPhi",&caloMetPhi,"CaloMETPhi");

 outTree->Branch("PFMET",&pfMet,"PFMET");
 outTree->Branch("PFMETPhi",&pfMetPhi,"PFMETPhi");

 outTree->Branch("TCMET",&tcMet,"TCMET");
 outTree->Branch("TCMETPhi",&tcMetPhi,"TCMETPhi");

 //FOR MONTE CARLO:
  if(!_isData){
    //generator level information
    outTree->Branch("procID",&procID,"procID/I");
    outTree->Branch("qScale",&qScale,"qScale/F");
    outTree->Branch("nPU",&nPu,"nPU/F");
    outTree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");

    ///gen electron, muon,photon
    outTree->Branch("nGenPho",&nGenPho,"nGenPho/I");
    outTree->Branch("GenPhotons",&GenPhotons);

    outTree->Branch("nGenEle",&nGenEle,"nGenEle/I");
    outTree->Branch("GenElectrons",&GenElectrons);
    
    outTree->Branch("nGenMu",&nGenMu,"nGenMu/I");
    outTree->Branch("GenMuons",&GenMuons);
    //generator level higgs
    outTree->Branch("nGenHiggs",&nGenHiggs,"nGenHiggs/I");
    outTree->Branch("GenHiggs",&GenHiggs);
    //other particles
    outTree->Branch("nGenOthers",&nGenOthers,"nGenOthers/I");
    outTree->Branch("GenOthers",&GenOthers);
  }

}
