#include "CommonTools/include/Utils.hh"
#include "include/SuperClustersEESelection.hh"
#include "include/JetCounter.hh"

#include "cajun/json/reader.h"
#include "cajun/json/elements.h"

#include <iostream>
#include <fstream>

using namespace std;

SuperClustersEESelection::SuperClustersEESelection(TTree *tree)
  : Vecbos(tree) {

  // default do not check on mc-truth
  m_signal = all;

  // common kinematic selections
  std::string theConfigDir = "config/superclustersel/";
  std::string fileCuts     = theConfigDir + "cuts.txt";
  std::string fileSwitches = theConfigDir + "switches.txt";

  _selection = new Selection(fileCuts,fileSwitches);
  _counters  = new Counters();
  ConfigSelection(_selection,_counters);

  // To read good run list!
  if (_selection->getSwitch("goodRunLS") && _selection->getSwitch("isData")) {
    std::string goodRunGiasoneFile = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

  WToENuDecay = 0;
  isData_ = _selection->getSwitch("isData");
  if(isData_) mcevent.SetData(true);
}

SuperClustersEESelection::~SuperClustersEESelection() {

  delete _counters;
  delete _selection;
}

void SuperClustersEESelection::ConfigSelection(Selection* selection, Counters* counters) {
  
  selection->addSwitch("isData");
  selection->addSwitch("mcTruth");
  selection->addSwitch("goodRunLS");
  selection->addSwitch("trigger");
  selection->addCut("eventRange");
  selection->addCut("etaSCAcc");
  selection->addCut("ptSCAcc");
  selection->addCut("spikeFraction");
  selection->addCut("seedTime");
  selection->addCut("HoE");
  selection->addCut("ECALisol");
  selection->addCut("etaCaloJetAcc");
  selection->addCut("etCaloJetAcc");
  selection->addCut("caloJetConeWidth");
  selection->addCut("etaPFJetAcc");
  selection->addCut("etPFJetAcc");
  selection->addCut("PFJetConeWidth");
  selection->addCut("Mt");
  selection->summary();

  counters->SetTitle("EVENT COUNTER");
  counters->AddVar("event");
  counters->AddVar("mcTruth");
  counters->AddVar("trigger");
  counters->AddVar("recoSCs");
  counters->AddVar("accSCs");
  counters->AddVar("spikeFraction");
  counters->AddVar("seedTime");
  counters->AddVar("idSCs");
  counters->AddVar("ECALisol");
  counters->AddVar("Mt");
  counters->AddVar("recoGSF");
  counters->AddVar("fullSelection");
}

void SuperClustersEESelection::Loop() {

  if(fChain == 0) return;

  char namefile[500];
  sprintf(namefile,"%s-SCStudy.root",_prefix);
  fileOut_ = TFile::Open(namefile, "recreate");
  this->createOutTree();

  cout << "Requested signal type = " << m_signal << endl;

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Total number of entries in the chain = " << nentries << std::endl;
  int maxEvents = nentries;
  if(_selection->getSwitch("eventRange")) {
    maxEvents = (int) _selection->getUpperCut("eventRange");
    cout << "WARNING! switch eventRange ON! Will run only on the first " << maxEvents << " events!" << endl;
  }
  
  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  for (Long64_t jentry=0; jentry<maxEvents;jentry++) {
    
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    reloadTriggerMask();

    // Good Run selection
    if (isData_ && _selection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun = runNumber;
        lastLumi = lumiBlock;
        std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (isData_ && _selection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    _counters->IncrVar("event");

    int indexMcEleWToENu = -1;
    if ( !isData_ ) { 
      mcevent.LoadDecay(nMc,idMc,mothMc);
      mcevent.LoadMomentum(pMc,energyMc,thetaMc,phiMc);
      // check if the decay is a W->enu prompt
      indexMcEleWToENu = mcevent.indexEleWPrompt();
      WToENuDecay = (indexMcEleWToENu >-1 ) ? 1 : 0;
    }

    // MC truth
    int mctruth = 0;
    if( m_signal == wjets ) mctruth = WToENuDecay;
    else if( m_signal == wother ) mctruth = !WToENuDecay;
    else mctruth = 1;
    if( _selection->getSwitch("mcTruth") && !mctruth ) continue;
    _counters->IncrVar("mcTruth");

    // trigger 
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();

    if( _selection->getSwitch("trigger") && !passedHLT) continue;
    _counters->IncrVar("trigger");

    if( nSC==0 ) continue;
    _counters->IncrVar("recoSCs");

    // loop over reconstructed ECAL superclusters and choose the accepted ones
    std::vector<int> acceptSCs;
    for(int isc=0; isc<nSC; isc++) {
      TVector3 pClu;
      pClu.SetMagThetaPhi(energySC[isc],thetaSC[isc],phiSC[isc]);
      if( _selection->getSwitch("etaSCAcc") && !_selection->passCut("etaSCAcc",fabs(etaSC[isc])) ) continue;
      bool isInEcalFiducial = false;
      if ( (fabs(etaSC[isc]) < 1.4442) || ( fabs(etaSC[isc])>1.560 && fabs(etaSC[isc])<2.5 ) ) isInEcalFiducial = true;
      if (!isInEcalFiducial) continue;
      if( _selection->getSwitch("ptSCAcc") && !_selection->passCut("ptSCAcc",pClu.Pt()) ) continue;
      acceptSCs.push_back(isc);
    }
    if (acceptSCs.size()>0) _counters->IncrVar("accSCs");

    // remove spikes... swiss-cross
    std::vector<int> noSpikes1;
    for(int i=0; i<(int)acceptSCs.size(); i++) {
      int isc  = acceptSCs[i];
      float e1 = eMaxSC[isc];
      float e4SwissCross = e4SwissCrossSC[isc];
      if(_selection->getSwitch("spikeFraction") && !_selection->passCut("spikeFraction", 1.0-e4SwissCross/e1) ) continue;
      noSpikes1.push_back(isc); 
    }
    if(noSpikes1.size()>0) _counters->IncrVar("spikeFraction");
    
    // remove spikes... timing
    std::vector<int> noSpikes;
    for(int i=0; i<(int)noSpikes1.size(); i++) {
      int isc = noSpikes1[i];
      float theSeedTime = timeSC[isc]; 
      if(_selection->getSwitch("seedTime") && !_selection->passCut("seedTime", theSeedTime) ) continue;  
      noSpikes.push_back(isc); 
    }

    // loop over accepted SCs not potential spikes and choose the hardest one
    int bestSC=-1;
    float maxEt=-1.0;
    for(int i=0; i<(int)noSpikes.size(); i++) {
      int isc = noSpikes[i];
      TVector3 pClu;
      pClu.SetMagThetaPhi(energySC[isc],thetaSC[isc],phiSC[isc]);
      if(pClu.Pt()>maxEt) {
        bestSC = isc;
        maxEt=pClu.Pt();
      }
    }
    if( bestSC<0 ) continue;
    _counters->IncrVar("seedTime");
    
    TVector3 pClu;
    pClu.SetMagThetaPhi(energySC[bestSC],thetaSC[bestSC],phiSC[bestSC]);

    std::vector<TVector3> pCluT3V;
    pCluT3V.push_back( pClu );

    // H/E cut fix: H/E in ntuples in H....
    float theHoE = hOverESC[bestSC]/energySC[bestSC];
    if( _selection->getSwitch("HoE") && !_selection->passCut("HoE",theHoE) ) continue;
    _counters->IncrVar("idSCs");

    // ECAL based isolation cut  [ isolation taken from gamma: if the corresponding photon is not recoed, keep the event anyway ]
    if (ecalRecHitSumEtConeDR04SC[bestSC]>-800) {
      if( _selection->getSwitch("ECALisol") && !_selection->passCut("ECALisol",ecalRecHitSumEtConeDR04SC[bestSC]) ) continue;
    }
    _counters->IncrVar("ECALisol");

    // count calojets for jet veto
    std::vector<Jet> _theAK5CaloJets = GetCorrJets(); 
    JetCounter goodCaloJets_counter(_theAK5CaloJets); 
    goodCaloJets_counter.SetThresholds(_selection->getLowerCut("etCaloJetAcc"), _selection->getUpperCut("etaCaloJetAcc"));
    goodCaloJets_counter.SetDistance(_selection->getUpperCut("caloJetConeWidth"));
    goodCaloJets_counter.SetParticlesToRemove(pCluT3V);
    std::vector<Jet> goodCaloJets_;
    goodCaloJets_ = goodCaloJets_counter.getGoodJets();
    int howManyCaloJets = goodCaloJets_.size();

    // count PFjets for jet veto
    std::vector<Jet> _theAK5PFJets = GetCorrPFJets();
    JetCounter goodPFJets_counter(_theAK5PFJets); 
    goodPFJets_counter.SetThresholds(_selection->getLowerCut("etPFJetAcc"), _selection->getUpperCut("etaPFJetAcc"));
    goodPFJets_counter.SetDistance(_selection->getUpperCut("PFJetConeWidth"));
    goodPFJets_counter.SetParticlesToRemove(pCluT3V);
    std::vector<Jet> goodPFJets_;
    goodPFJets_ = goodPFJets_counter.getGoodJets();
    int howManyPFJets = goodPFJets_.size();

    // transverse mass and met
    TVector3 p3PFMet(pxPFMet[0],pyPFMet[0],0.0);
    TVector3 pT3SC(pClu.Px(),pClu.Py(),0.0);
    float WPFmT = sqrt(2 * pClu.Pt() * p3PFMet.Mag() * (1-cos(pT3SC.Angle(p3PFMet))) );   
    float WPFmeT = GetPt(pxPFMet[0],pyPFMet[0]);

    if( _selection->getSwitch("Mt") && !_selection->passCut("Mt",WPFmT) ) continue;
    _counters->IncrVar("Mt");

    _counters->IncrVar("fullSelection");
    
    // counting the number of SC above threshold and in the acceptance
    int nSCok = 0;
    for (int iSC=0; iSC<nSC; iSC++) {
      float etSC = energySC[iSC]*sin(thetaSC[iSC]);
      if( !_selection->passCut("ptSCAcc",etSC) ) continue;	  
      if( !_selection->passCut("etaSCAcc",fabs(etaSC[iSC])) ) continue;
      bool isInEcalFiducial = false;
      if( (fabs(etaSC[iSC]) < 1.4442) || ( fabs(etaSC[iSC])>1.560 && fabs(etaSC[iSC])<2.5 ) ) isInEcalFiducial = true;
      if( !isInEcalFiducial ) continue;
      nSCok++;
    }

    if (bestSC<=-1) { cout << "chiara: problema enorme" << endl; }

    // look for the GSF electron
    int matchedSC=0;
    for(int iele=0; iele<nEle; iele++) {
      int eleSc = superClusterIndexEle[iele];
      if( eleSc == bestSC ) {
        matchedSC = 1;
        break;
      }
    }
    if( matchedSC ) _counters->IncrVar("recoGSF");

    // look for a KF track
    int theSelKFTrack     = -1;
    float deltaMinKFTrack = 999.;
    if (bestSC>-1) {
      for (int theTrack=0; theTrack<nTrack; theTrack++) {
	TVector3 pTrack(pxTrack[theTrack],pyTrack[theTrack],pzTrack[theTrack]);
	float thisDeltaR = pClu.DeltaR(pTrack);
	if (fabs(thisDeltaR)<fabs(deltaMinKFTrack)) {
	  deltaMinKFTrack = fabs(thisDeltaR);
	  theSelKFTrack   = theTrack;
	}
      }
    }

    // look for a GSF track 
    int theSelGsfTrack     = -1;
    float deltaMinGsfTrack = 999.;
    if (bestSC>-1) {
      for (int theGsfTrack=0; theGsfTrack<nGsfTrack; theGsfTrack++) {
	TVector3 pGsfTrack(pxGsfTrack[theGsfTrack],pyGsfTrack[theGsfTrack],pzGsfTrack[theGsfTrack]);
	float thisDeltaR = pClu.DeltaR(pGsfTrack);
	if (fabs(thisDeltaR)<fabs(deltaMinGsfTrack)) {
	  deltaMinGsfTrack = fabs(thisDeltaR);
	  theSelGsfTrack   = theGsfTrack;
	}
      }
    }

    // look for a closer PF electron
    int theSelPFele = -1;
    float deltaMin  = 999.;
    if (bestSC>-1) {
      for (int thePFEle=0; thePFEle<nPFEle; thePFEle++) {
	TVector3 pPFelectron(pxPFEle[thePFEle],pyPFEle[thePFEle],pzPFEle[thePFEle]);
	float thisDeltaR = pClu.DeltaR(pPFelectron);
	if (fabs(thisDeltaR)<fabs(deltaMin)) {
	  deltaMin    = fabs(thisDeltaR);
	  theSelPFele = thePFEle;
	}
      }
    }
  
    // fill the output tree
    myPfElePt     = -5000.;
    myPfEleEta    = -5000.;
    myPfElePhi    = -5000.;
    myPfEleDeltaR = -5000.;
    myPfEleMva    = -5000.;
    myNPFEle      = -5000;

    if (theSelPFele>=0) {
      TVector3 selPFelectron(pxPFEle[theSelPFele],pyPFEle[theSelPFele],pzPFEle[theSelPFele]);
      myPfElePt     = selPFelectron.Pt();
      myPfEleEta    = selPFelectron.Eta();
      myPfElePhi    = selPFelectron.Phi();
      myPfEleDeltaR = deltaMin;
      myPfEleMva    = MvaOutputPFEle[theSelPFele];
      myNPFEle      = nPFEle;
    }
    if (theSelPFele<0 && bestSC>-1) {
      myPfElePt     = -999.;
      myPfEleEta    = -999.;
      myPfElePhi    = -999.;
      myPfEleDeltaR = -999.;
      myPfEleMva    = -999.;
      myNPFEle      = nPFEle;
    }
    if (theSelPFele<0 && bestSC==-1) {
      myPfElePt     = -1.;
      myPfEleEta    = -1.;
      myPfElePhi    = -1.;
      myPfEleDeltaR = -1.;
      myPfEleMva    = -1.;
      myNPFEle      = nPFEle;
    }
    
    myMt  = WPFmT;
    myMet = WPFmeT;
    myEta = pClu.Eta();
    myPt  = pClu.Pt();
    myPhi = pClu.Phi();      
    
    myRecoGSF   = matchedSC;
    myNSCall    = nSC;
    myNSCok     = nSCok;
    myNCaloJets = howManyCaloJets;
    myNPFJets   = howManyPFJets;

    myHoE = theHoE;
    myE9  = e3x3SC[bestSC];
    myE25 = e5x5SC[bestSC];
    myCovIeIe = covIEtaIEtaSC[bestSC];
    myScBasedEcalSum03SC        = scBasedEcalSum03SC[bestSC]; 
    myEcalRecHitSumEtConeDR03SC = ecalRecHitSumEtConeDR03SC[bestSC];
    myHcalTowerSumEtConeDR03SC  = hcalTowerSumEtConeDR03SC[bestSC];
    myTrkSumPtSolidConeDR03SC   = trkSumPtSolidConeDR03SC[bestSC];

    myRunNumber   = runNumber;
    myEventNumber = eventNumber;

    myMatchedEG = 0.;
    if ( matchedSC ) myMatchedEG=1;
    if (!matchedSC ) myMatchedEG=-1;

    myMatchedPF = 0.;
    if (myPfEleDeltaR<0.15 && myPfEleDeltaR>-400) myMatchedPF = 1;
    if (myPfEleDeltaR>0.15 || myPfEleDeltaR<-400) myMatchedPF =-1;

    myMatchedKF = 0.;
    if (deltaMinKFTrack<0.15) myMatchedKF = 1;
    if (deltaMinKFTrack>0.15) myMatchedKF = -1;

    myMatchedGsf = 0.;
    if (deltaMinGsfTrack<0.15) myMatchedGsf = 1;
    if (deltaMinGsfTrack>0.15) myMatchedGsf = -1;

    outTree_->Fill();
  }

  outTree_->Write();
  fileOut_->Close();

}

void SuperClustersEESelection::setJsonGoodRunList(const string& jsonFilePath) {

  jsonFile=jsonFilePath;
}

void SuperClustersEESelection::fillRunLSMap() {
  
  if (jsonFile == "") {
    std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
    return;
  }

  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open()) {
    std::cout << "Unable to open file " << jsonFile << std::endl;
    return;
  }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);
  
  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun) {

    const json::Array& lsSegment = (*itRun).element;
    LSSegments thisRunSegments; 
    for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator) {

      json::Array lsSegment=(*lsIterator);
      json::Number lsStart=lsSegment[0];	   
      json::Number lsEnd=lsSegment[1];
      aLSSegment thisSegment;
      thisSegment.first=lsStart.Value();
      thisSegment.second=lsEnd.Value();
      thisRunSegments.push_back(thisSegment);
    }
    goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
  }
  

  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR) {
    std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
    for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
      std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] "; 
    std::cout << std::endl;
  }
}

bool SuperClustersEESelection::isGoodRunLS() {
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg) {
    if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
      return true;
  }
  return false;
}

void SuperClustersEESelection::createOutTree() {

  outTree_ = new TTree("T1","supercluser tree");
  outTree_->Branch("eta", &myEta, "eta/F");
  outTree_->Branch("pt",  &myPt,  "pt/F");
  outTree_->Branch("mt",  &myMt,  "mt/F");
  outTree_->Branch("met", &myMet, "met/F");
  outTree_->Branch("hoe", &myHoE, "hoe/F");
  outTree_->Branch("recoGSF",     &myRecoGSF,     "recoGSF/I");
  outTree_->Branch("pfElePt",     &myPfElePt,     "pfElePt/F");
  outTree_->Branch("pfEleEta",    &myPfEleEta,    "pfEleEta/F");
  outTree_->Branch("pfElePhi",    &myPfElePhi,    "pfElePhi/F");
  outTree_->Branch("pfEleDeltaR", &myPfEleDeltaR, "pfEleDeltaR/F");
  outTree_->Branch("pfEleMva",    &myPfEleMva,    "pfEleMva/F");
  outTree_->Branch("nPFEle",      &myNPFEle,      "nPFEle/I");
  outTree_->Branch("nSCall",      &myNSCall,      "nSCall/I");
  outTree_->Branch("nSCok",       &myNSCok,       "nSCok/I");
  outTree_->Branch("matchedEG",   &myMatchedEG,   "matchedEG/I");
  outTree_->Branch("matchedPF",   &myMatchedPF,   "matchedPF/I");
  outTree_->Branch("matchedKF",   &myMatchedKF,   "matchedKF/I");
  outTree_->Branch("matchedGsf",  &myMatchedGsf,  "matchedGsf/I");
  outTree_->Branch("nCaloJets",   &myNCaloJets,   "nCaloJets/I");
  outTree_->Branch("nPFJets",     &myNPFJets,     "nPFJets/I");
  outTree_->Branch("runNumber",   &myRunNumber,   "runNumber/I");
  outTree_->Branch("eventNumber", &myEventNumber, "eventNumber/I");
  outTree_->Branch("hoe",         &myHoE,         "hoe/F");  
  outTree_->Branch("phi",         &myPhi,         "phi/F");
  outTree_->Branch("E9",          &myE9,          "E9/F"); 
  outTree_->Branch("E25",         &myE25,         "E25/F"); 
  outTree_->Branch("covIeIe",     &myCovIeIe,     "covIeIe/F"); 
  outTree_->Branch("scBasedEcalSum03SC",        &myScBasedEcalSum03SC,        "scBasedEcalSum03SC/F");
  outTree_->Branch("ecalRecHitSumEtConeDR03SC", &myEcalRecHitSumEtConeDR03SC, "ecalRecHitSumEtConeDR03SC/F");
  outTree_->Branch("hcalTowerSumEtConeDR03SC",  &myHcalTowerSumEtConeDR03SC,  "hcalTowerSumEtConeDR03SC/F");
  outTree_->Branch("trkSumPtSolidConeDR03SC",   &myTrkSumPtSolidConeDR03SC,   "trkSumPtSolidConeDR03SC/F");
}

void SuperClustersEESelection::displayEfficiencies() {
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "+++ DETAILED EFFICIENCY +++ " << std::endl;
  
  char namefile[500];
  sprintf(namefile,"%s-SC-Counters.root",_prefix);
  
  _counters->Draw();
  _counters->Draw("mcTruth","event");
  _counters->Draw("trigger","mcTruth");
  _counters->Draw("recoSCs","trigger");
  _counters->Draw("accSCs","recoSCs");
  _counters->Draw("idSCs","accSCs");
  _counters->Draw("ECALisol","idSCs");
  _counters->Draw("Mt","ECALisol");
  _counters->Draw("fullSelection","mcTruth");
  _counters->Draw("recoGSF","fullSelection");
  
  _counters->Save(namefile,"recreate");

}
