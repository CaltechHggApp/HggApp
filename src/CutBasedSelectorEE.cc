#include "include/CutBasedSelectorEE.hh"
#include <iostream>


SelectorData::SelectorData() {
  signal = 0;
  inacceptance = 0;
  njetsMC = 0;
  njets = 0;
  passedHLT = false;
  matchedHLT = false;
  nRecoEle = 0;
  nAccEle = 0;
  nIdTightEle = 0;
  nIdLooseEle = 0;
  nIsolTightEle = 0;
  nIsolLooseEle = 0;
  nConvRejTightEle = 0;
  nConvRejLooseEle = 0;
  nPV = 0;
  nDzVertexEle = 0;
  nDxyVertexEle = 0;
  nMuons = 0;
  btagEVT = 0;
  mInv = 0;
  mhtMET = 0;
  met = 0;
  mt = 0;
  chargeMajorityMethod = false;
  foundAnyZ = false;
  mhtJet.reserve(6);
  for(int j=0; j<6; j++) mhtJet[j] = 0;
  weight = 0;
}

SelectorData::SelectorData(SelectorData &data) {
  signal = data.signal;
  inacceptance = data.inacceptance;
  njetsMC = data.njetsMC;
  njets = data.njets;
  passedHLT = data.passedHLT;
  matchedHLT = data.matchedHLT;
  nRecoEle = data.nRecoEle;
  nAccEle = data.nAccEle;
  nIdTightEle = data.nIdTightEle;
  nIdLooseEle = data.nIdLooseEle;
  nIsolTightEle = data.nIsolTightEle;
  nIsolLooseEle = data.nIsolLooseEle;
  nConvRejTightEle = data.nConvRejTightEle;
  nConvRejLooseEle = data.nConvRejLooseEle;
  nPV = data.nPV;
  nDzVertexEle = data.nDzVertexEle;
  nDxyVertexEle = data.nDxyVertexEle;
  nMuons = data.nMuons;
  btagEVT = data.btagEVT;
  mInv = data.mInv;
  mhtMET = data.mhtMET;
  chargeMajorityMethod = data.chargeMajorityMethod;
  met = data.met;
  mt = data.mt;
  foundAnyZ = data.foundAnyZ;
  mhtJet = data.mhtJet;
  weight = data.weight;
}

CutBasedSelectorEE::CutBasedSelectorEE(SelectorData data) {
  data_=data;
  mc_=0;
}

bool CutBasedSelectorEE::outputZ(Selection* commonSelection, Selection* selection, 
                                 Counters* inclusiveCounter, 
                                 Counters (*mcJetBinCounter)[6], Counters (*recoJetBinCounter)[6]) {
  // all events
  if( inclusiveCounter ) inclusiveCounter->IncrVar("event",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("event",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("event",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("event",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("event",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("event",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("event",data_.weight);
  }
  
  
  // MC truth
  if(!commonSelection->getSwitch("isData") && !data_.signal ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("mcTruth",data_.weight);
  
  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("mcTruth",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("mcTruth",data_.weight);
  } else {    
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("mcTruth",data_.weight);
    else (*mcJetBinCounter)[4].IncrVar("mcTruth",data_.weight);
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("mcTruth",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("mcTruth",data_.weight);
  }


  // within MC acceptance  
  if(!commonSelection->getSwitch("isData") && commonSelection->getSwitch("inacceptance") && !data_.inacceptance ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("inacceptance",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("inacceptance",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("inacceptance",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("inacceptance",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("inacceptance",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("inacceptance",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("inacceptance",data_.weight);
  }


  // hlt
  if(selection->getSwitch("trigger") && !data_.passedHLT ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("trigger",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("trigger",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("trigger",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("trigger",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("trigger",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("trigger",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("trigger",data_.weight);
  }  


  // electrons multiplicity
  if (selection->getSwitch("nRecoEles") && !selection->passCut("nRecoEles", data_.nRecoEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nRecoEles",data_.weight);
  
  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nRecoEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nRecoEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nRecoEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nRecoEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nRecoEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nRecoEles",data_.weight);
  }

  if (selection->getSwitch("nAccEles") && !selection->passCut("nAccEles", data_.nAccEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nAccEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nAccEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nAccEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nAccEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nAccEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nAccEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nAccEles",data_.weight);
  }


  bool passEleId=false;
  if(!commonSelection->getSwitch("asymmetricElectrons") ) {
    if(data_.nIdTightEle>=2) passEleId=true;  
  } else {
    if(data_.nIdTightEle>=2) passEleId=true;
    else if(data_.nIdTightEle==1 && data_.nIdLooseEle>=1) passEleId=true;
  }

  if (selection->getSwitch("nIdEles") && !passEleId) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nIdEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nIdEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nIdEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nIdEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nIdEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nIdEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nIdEles",data_.weight);
  }

  bool passIsol=false;
  if(!commonSelection->getSwitch("asymmetricElectrons") ) {
    if(data_.nIsolTightEle>=2) passIsol=true;  
  } else {
    if(data_.nIsolTightEle>=2) passIsol=true;
    else if(data_.nIsolTightEle==1 && data_.nIsolLooseEle>=1) passIsol=true;
  }
  
  if (selection->getSwitch("nIsolEles") && !passIsol) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nIsolEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nIsolEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nIsolEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nIsolEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nIsolEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nIsolEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nIsolEles",data_.weight);
  }

  bool passConvRej=false;
  if(!commonSelection->getSwitch("asymmetricElectrons") ) {
    if(data_.nConvRejTightEle>=2) passConvRej=true;  
  } else {
    if(data_.nConvRejTightEle>=2) passConvRej=true;
    else if(data_.nConvRejTightEle==1 && data_.nConvRejLooseEle>=1) passConvRej=true;
  }
  
  if (selection->getSwitch("nConvRejEles") && !passConvRej) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nConvRejEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nConvRejEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nConvRejEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nConvRejEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nConvRejEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nConvRejEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nConvRejEles",data_.weight);
  }


  if (selection->getSwitch("nPrimaryVertices") && !selection->passCut("nPrimaryVertices", data_.nPV)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nPrimaryVertices",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nPrimaryVertices",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nPrimaryVertices",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nPrimaryVertices",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nPrimaryVertices",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nPrimaryVertices",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nPrimaryVertices",data_.weight);
  }


  if (selection->getSwitch("nDzVertexEles") && !selection->passCut("nDzVertexEles", data_.nDzVertexEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nDzVertexEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nDzVertexEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nDzVertexEles",data_.weight);  
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nDzVertexEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nDzVertexEles",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nDzVertexEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nDzVertexEles",data_.weight);
  }


  if (selection->getSwitch("nDxyVertexEles") && !selection->passCut("nDxyVertexEles", data_.nDxyVertexEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nDxyVertexEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nDxyVertexEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nDxyVertexEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nDxyVertexEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nDxyVertexEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nDxyVertexEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nDxyVertexEles",data_.weight);
  }

  if(selection->getSwitch("matchedHLT") && !data_.matchedHLT) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("matchedHLT",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("matchedHLT",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("matchedHLT",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("matchedHLT",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("matchedHLT",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("matchedHLT",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("matchedHLT",data_.weight);
  }


  // additional muon veto (from ttbar decays)
  if (selection->getSwitch("nMuons") && !selection->passCut("nMuons", data_.nMuons)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nMuons",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nMuons",data_.weight);
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nMuons",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nMuons",data_.weight);
    else (*mcJetBinCounter)[4].IncrVar("nMuons",data_.weight);
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nMuons",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nMuons",data_.weight);
  }


  // B veto
  if (selection->getSwitch("btagEVT") && !selection->passCut("btagEVT", data_.btagEVT) ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("btagEVT",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("btagEVT",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("btagEVT",data_.weight);
  } else{ 
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("btagEVT",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("btagEVT",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("btagEVT",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("btagEVT",data_.weight);
  }


  // ee invariant mass
  if (selection->getSwitch("meeCut") && !selection->passCut("meeCut", data_.mInv) ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("meeCut",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("meeCut",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("meeCut",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("meeCut",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("meeCut",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("meeCut",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("meeCut",data_.weight);
  }


  // event shape variable
  if (selection->getSwitch("MHTphiMET") && !selection->passCut("MHTphiMET", data_.mhtMET) ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("MHTphiMET",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("MHTphiMET",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("MHTphiMET",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("MHTphiMET",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("MHTphiMET",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("MHTphiMET",data_.weight);
    else(*recoJetBinCounter)[4].IncrVar("MHTphiMET",data_.weight); 
  }


  if( inclusiveCounter ) inclusiveCounter->IncrVar("fullSelection",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("fullZeeSelection",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("fullZeeSelection",data_.weight);
  } else {

    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("fullZeeSelection",data_.weight); 
    else(*mcJetBinCounter)[4].IncrVar("fullZeeSelection",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("fullZeeSelection",data_.weight);
    else(*recoJetBinCounter)[4].IncrVar("fullZeeSelection",data_.weight);
  }

  return true;
}

bool CutBasedSelectorEE::outputW(Selection* commonSelection, Selection* selection, 
                                 Counters* inclusiveCounter, 
                                 Counters (*mcJetBinCounter)[6], Counters (*recoJetBinCounter)[6]) {

  // all events
  if( inclusiveCounter ) inclusiveCounter->IncrVar("event",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("event",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("event",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("event",data_.weight); 
    else(*mcJetBinCounter)[4].IncrVar("event",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("event",data_.weight);
    else(*recoJetBinCounter)[4].IncrVar("event",data_.weight);
  }  


  // MC truth
  if(!commonSelection->getSwitch("isData") && !data_.signal ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("mcTruth",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("mcTruth",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("mcTruth",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("mcTruth",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("mcTruth",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("mcTruth",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("mcTruth",data_.weight);
  }  


  // within MC acceptance  
  if(!commonSelection->getSwitch("isData") && commonSelection->getSwitch("inacceptance") && !data_.inacceptance ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("inacceptance",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("inacceptance",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("inacceptance",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("inacceptance",data_.weight); 
    else(*mcJetBinCounter)[4].IncrVar("inacceptance",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("inacceptance",data_.weight);
    else(*recoJetBinCounter)[4].IncrVar("inacceptance",data_.weight);
  }


  // hlt
  if(selection->getSwitch("trigger") && !data_.passedHLT ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("trigger",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("trigger",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("trigger",data_.weight);
  } else {

    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("trigger",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("trigger",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("trigger",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("trigger",data_.weight);
  }


  // electrons multiplicity
  if (selection->getSwitch("nRecoEles") && !selection->passCut("nRecoEles", data_.nRecoEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nRecoEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nRecoEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nRecoEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("nRecoEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nRecoEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nRecoEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nRecoEles",data_.weight);
  }

  if (selection->getSwitch("nAccEles") && !selection->passCut("nAccEles", data_.nAccEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nAccEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nAccEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nAccEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("nAccEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nAccEles",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nAccEles",data_.weight);
    else(*recoJetBinCounter)[4].IncrVar("nAccEles",data_.weight);
  }

  if (selection->getSwitch("nIdEles") && !selection->passCut("nIdEles", data_.nIdTightEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nIdEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nIdEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nIdEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("nIdEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nIdEles",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nIdEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nIdEles",data_.weight);
  }

  if (selection->getSwitch("nIsolEles") && !selection->passCut("nIsolEles", data_.nIsolTightEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nIsolEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nIsolEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nIsolEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("nIsolEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nIsolEles",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nIsolEles",data_.weight);
    else(*recoJetBinCounter)[4].IncrVar("nIsolEles",data_.weight);
  }

  if (selection->getSwitch("nConvRejEles") && !selection->passCut("nConvRejEles", data_.nConvRejTightEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nConvRejEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nConvRejEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nConvRejEles",data_.weight);
  } else { 
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("nConvRejEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nConvRejEles",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nConvRejEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nConvRejEles",data_.weight);
  }

  if(selection->getSwitch("ZVeto") && data_.foundAnyZ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("ZVeto",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("ZVeto",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("ZVeto",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("ZVeto",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("ZVeto",data_.weight);  
    if ( data_.njets<4 ) (*recoJetBinCounter)[data_.njets].IncrVar("ZVeto",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("ZVeto",data_.weight);
  }

  if (selection->getSwitch("nPrimaryVertices") && !selection->passCut("nPrimaryVertices", data_.nPV)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nPrimaryVertices",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nPrimaryVertices",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nPrimaryVertices",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("nPrimaryVertices",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nPrimaryVertices",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nPrimaryVertices",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nPrimaryVertices",data_.weight);
  }

  if (selection->getSwitch("nDzVertexEles") && !selection->passCut("nDzVertexEles", data_.nDzVertexEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nDzVertexEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nDzVertexEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nDzVertexEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("nDzVertexEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nDzVertexEles",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nDzVertexEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nDzVertexEles",data_.weight);
  }  

  if (selection->getSwitch("nDxyVertexEles") && !selection->passCut("nDxyVertexEles", data_.nDxyVertexEle)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nDxyVertexEles",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nDxyVertexEles",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nDxyVertexEles",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nDxyVertexEles",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("nDxyVertexEles",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nDxyVertexEles",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nDxyVertexEles",data_.weight);
  }

  // charge
  if(selection->getSwitch("charge") && !data_.chargeMajorityMethod ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("charge",data_.weight);
  
  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("charge",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("charge",data_.weight);
  } else {
    
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("charge",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("charge",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("charge",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("charge",data_.weight);
  }

  if(selection->getSwitch("matchedHLT") && !data_.matchedHLT) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("matchedHLT",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("matchedHLT",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("matchedHLT",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("matchedHLT",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("matchedHLT",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("matchedHLT",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("matchedHLT",data_.weight);
  }


  // additional muon veto (from ttbar decays)
  if (selection->getSwitch("nMuons") && !selection->passCut("nMuons", data_.nMuons)) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("nMuons",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("nMuons",data_.weight);
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("nMuons",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("nMuons",data_.weight);
    else (*mcJetBinCounter)[4].IncrVar("nMuons",data_.weight); 
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("nMuons",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("nMuons",data_.weight);
  }


  // B veto
  if (selection->getSwitch("btagEVT") && !selection->passCut("btagEVT", data_.btagEVT) ) return false;
  if( inclusiveCounter ) inclusiveCounter->IncrVar("btagEVT",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("btagEVT",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("btagEVT",data_.weight);
  } else {
    if ( data_.njetsMC<4 ) (*mcJetBinCounter)[data_.njetsMC].IncrVar("btagEVT",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("btagEVT",data_.weight);  
    if ( data_.njets<4) (*recoJetBinCounter)[data_.njets].IncrVar("btagEVT",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("btagEVT",data_.weight);
  }


  // met value
  if (selection->getSwitch("metCut") && !selection->passCut("metCut", data_.met) ) return false;  
  if( inclusiveCounter ) inclusiveCounter->IncrVar("metCut",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("metCut",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("metCut",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("metCut",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("metCut",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("metCut",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("metCut",data_.weight);
  }


  // transverse mass value
  if (selection->getSwitch("transvMassCut") && !selection->passCut("transvMassCut", data_.mt) ) return false;  
  if( inclusiveCounter ) inclusiveCounter->IncrVar("transvMassCut",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("transvMassCut",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("transvMassCut",data_.weight);
  } else{ 
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("transvMassCut",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("transvMassCut",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("transvMassCut",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("transvMassCut",data_.weight);
  }


  // event shape variable:
  // the variable is different depending on the jet multiplicity. So the selection cannot be done at 
  // this step because we dump inclusive trees
  // The real cut has to be done when producing 1 dataset for each specific multiplicity
  // For the efficiency, it can be done here because we have a counter for each jet multiplicity
  if( inclusiveCounter ) inclusiveCounter->IncrVar("MHTphiJet",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) { 
      if( selection->getSwitch("MHTphiJet") && selection->passCut("MHTphiJet", data_.mhtJet[jet]) ) 
	(*mcJetBinCounter)[jet].IncrVar("MHTphiJet",data_.weight); 
    }
    for(int jet=0; jet<=data_.njets && jet<6; jet++) {
      if( selection->getSwitch("MHTphiJet") && selection->passCut("MHTphiJet", data_.mhtJet[jet]) )
	(*recoJetBinCounter)[jet].IncrVar("MHTphiJet",data_.weight);
    }
  } else { // this is not correct! In any case it is off now....
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("MHTphiJet",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("MHTphiJet",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("MHTphiJet",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("MHTphiJet",data_.weight);
  }


  if( inclusiveCounter ) inclusiveCounter->IncrVar("fullSelection",data_.weight);

  if(commonSelection->getSwitch("inclusiveCounting")) {
    for(int jet=0; jet<=data_.njetsMC && jet<6; jet++) (*mcJetBinCounter)[jet].IncrVar("fullSelection",data_.weight); 
    for(int jet=0; jet<=data_.njets && jet<6; jet++) (*recoJetBinCounter)[jet].IncrVar("fullSelection",data_.weight);
  } else {
    if ( data_.njetsMC<4 )(*mcJetBinCounter)[data_.njetsMC].IncrVar("fullSelection",data_.weight); 
    else (*mcJetBinCounter)[4].IncrVar("fullSelection",data_.weight);  
    if ( data_.njets<4 )(*recoJetBinCounter)[data_.njets].IncrVar("fullSelection",data_.weight);
    else (*recoJetBinCounter)[4].IncrVar("fullSelection",data_.weight);
  }

  return true;
}
