#include <VecbosEGObject.hh>
#include "CommonTools/include/Utils.hh"

#define debugEGObject 0
#include <iostream>
using namespace std;

VecbosBC::VecbosBC(){

}

VecbosBC::VecbosBC(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosBC::Init(VecbosBase* o, int i){
  if(i>o->nBC){
    index = -1;
    return;
  }
  index   = i;
  energy  = o->energyBC[i];
  eta     = o->etaBC[i];
  phi     = o->phiBC[i];

  
  etaCrystal = o->etaCrystalBC[i];
  phiCrystal = o->phiCrystalBC[i];
  iEta = o->iEtaBC[i];
  iPhi = o->iPhiBC[i];
  thetaTilt = o->thetaTiltBC[i];
  phiTilt   = o->phiTiltBC[i];
};

VecbosPFBC::VecbosPFBC(){

}

VecbosPFBC::VecbosPFBC(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosPFBC::Init(VecbosBase* o, int i){
  if(i>o->nPFBC) return;
  index   = i;
  energy  = o->energyPFBC[i];
  eta     = o->etaPFBC[i];
  phi     = o->phiPFBC[i];
  
  etaCrystal = o->etaCrystalPFBC[i];
  phiCrystal = o->phiCrystalPFBC[i];
  iEta = o->iEtaPFBC[i];
  iPhi = o->iPhiPFBC[i];
  thetaTilt = o->thetaTiltPFBC[i];
  phiTilt   = o->phiTiltPFBC[i];

};


struct sort_pred {
  bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right) {
    return left.second > right.second;
  }
};

VecbosSC::VecbosSC(){}

VecbosSC::VecbosSC(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosSC::Init(VecbosBase* o, int i){
  if(i>o->nSC || i<0){
    index = -1;
    return;
  }    
  index    = i;
  energy   = o->energySC[i];
  esEnergy = o->esEnergySC[i];
  eta      = o->etaSC[i];
  phi      = o->phiSC[i];
  e3x3     = o->e3x3SC[i];
  e5x5     = o->e5x5SC[i];

  e3x1     = o->e3x1SC[i];
  e1x3     = o->e1x3SC[i];
  e4x4     = o->e4x4SC[i];
  eMax     = o->eMaxSC[i];
  e2x2     = o->e2x2SC[i];
  e2nd     = o->e2ndSC[i];
  e1x5     = o->e1x5SC[i];
  e2x5Max  = o->e2x5MaxSC[i];
  e2x5Left = o->e2x5LeftSC[i];
  e2x5Right = o->e2x5RightSC[i];
  e2x5Top = o->e2x5TopSC[i];
  e2x5Bottom = o->e2x5BottomSC[i];
  
  eLeft    = o->eLeftSC[i];
  eRight   = o->eRightSC[i];
  eTop     = o->eTopSC[i];
  eBottom  = o->eBottomSC[i];


  sigmaIEtaIEta = sqrt(o->covIEtaIEtaSC[i]);
  sigmaIEtaIPhi = o->covIEtaIPhiSC[i];
  sigmaIPhiIPhi = sqrt(o->covIPhiIPhiSC[i]);

  esEffSigRR = TMath::Sqrt( TMath::Power(o->esEffsIxIxSC[i],2) +
			    TMath::Power(o->esEffsIyIySC[i],2) );

  CaloPos.SetXYZ(o->xPosSC[i],o->yPosSC[i],o->zPosSC[i]);

  rawE     = o->rawEnergySC[i];
  phiWidth = o->phiWidthSC[i];
  etaWidth = o->etaWidthSC[i];
  HoverE   = o->hOverESC[i];
  r9       = e3x3/rawE;
  s4ratio  = e2x2/e5x5;
  //get the basic clusters


  nBCs = o->nBCSC[i];
  std::vector< std::pair<int,float> > indexEnergyMap;
  float maxE=-1;
  int maxI=-1;
  for(int j=0; j<o->nBC; j++){
    if(o->indexSCBC[j] == i){ // the basic cluster points to this SC
      indexEnergyMap.push_back(std::pair<int,float>(j,o->energyBC[j]) );
      if(o->energyBC[j] > maxE){
	maxE = o->energyBC[j];
	maxI = j;
      }
    }
  }
  //std::sort(indexEnergyMap.begin(),indexEnergyMap.end(),sort_pred()); // sort the BCs by energy
  //for(int j=0;j<indexEnergyMap.size();j++){
  //  basicClusters.push_back(VecbosBC(o,indexEnergyMap.at(j).first));
  //}
  
  BCSeed.Init(o,maxI);
  BCSeed.energy = o->seedClusterEnergySC[i];
 }

/*
SCInfo VecbosSC::getStruct(){
  SCInfo out = {
  index,
  0,   // we set the BC pointer seperately
  energy,
  esEnergy, // in the trees, but not yet in VecbosBase
 eta,
  phi,
  e3x3,
  e5x5,
  sigmaIEtaIEta,
  sigmaIEtaIPhi,
  sigmaIPhiIPhi,
  rawE,
  phiWidth,
  etaWidth,
  HoverE,
  this->r9()
  };

  out.BCs[0]=-1.;   out.BCs[1]=-1.;   out.BCs[2]=-1.;   out.BCs[3]=-1.; //initialize
  if(basicClusters.size() >0) out.BCs[0] = basicClusters.at(0).index;
  if(basicClusters.size() >1) out.BCs[1] = basicClusters.at(1).index;
  if(basicClusters.size() >2) out.BCs[2] = basicClusters.at(basicClusters.size()-1).index;
  if(basicClusters.size() >0) out.BCs[0] = basicClusters.at(basicClusters.size()-2).index;

}
*/
VecbosPFSC::VecbosPFSC(){}

VecbosPFSC::VecbosPFSC(VecbosBase* o, int i){
  this->Init(o,i);
}
VecbosPFSC::VecbosPFSC(VecbosBase* o, float etaSC,float phiSC){
  this->Init(o,etaSC,phiSC);
}
void VecbosPFSC::Init(VecbosBase* o, float etaSC,float phiSC){
  const float maxDR = 0.3;

  float DR=9999;  int index = -1;
  for(int i=0;i<o->nPhoPFSC;i++){
    float thisDR = DeltaR(etaSC,o->etaPhoPFSC[i],phiSC,o->phiPhoPFSC[i]);
    if(thisDR < DR && thisDR < maxDR){
      DR = thisDR;
      index = i;      
    }
  }
  this->Init(o,index);

}
void VecbosPFSC::Init(VecbosBase* o, int i){
  if(i>o->nPhoPFSC || i<0){
    index = -1;
    return;
  }
  index    = i;
  energy   = o->energyPhoPFSC[i];
  esEnergy = o->esEnergyPhoPFSC[i];
  eta      = o->etaPhoPFSC[i];
  phi      = o->phiPhoPFSC[i];
  e3x3     = o->e3x3PhoPFSC[i];
  e5x5     = o->e5x5PhoPFSC[i];
  sigmaIEtaIEta = o->covIEtaIEtaPhoPFSC[i];
  sigmaIEtaIPhi = o->covIEtaIPhiPhoPFSC[i];
  sigmaIPhiIPhi = o->covIPhiIPhiPhoPFSC[i];

  CaloPos.SetXYZ(o->xPosPFSC[i],o->yPosPFSC[i],o->zPosPFSC[i]);


  rawE     = o->rawEnergyPhoPFSC[i];
  phiWidth = o->phiWidthPhoPFSC[i];
  etaWidth = o->etaWidthPhoPFSC[i];
  HoverE   = o->hOverEPhoPFSC[i];
  r9       = e3x3/rawE;

  //get the basic clusters
  std::vector< std::pair<int,float> > indexEnergyMap;
  for(int j=0; j<o->nPFBC; j++){
    if(o->indexSCPFBC[j] == i){ // the basic cluster points to this SC
      indexEnergyMap.push_back(std::pair<int,float>(j,o->energyPFBC[j]) );
    }
  }
  std::sort(indexEnergyMap.begin(),indexEnergyMap.end(),sort_pred()); // sort the BCs by energy
  for(int j=0;j<indexEnergyMap.size();j++){
    pfClusters.push_back(VecbosPFBC(o,indexEnergyMap.at(j).first));
  }

 }

 VecbosConversion::VecbosConversion(){}

VecbosConversion::VecbosConversion(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosConversion::Init(VecbosBase* o, int i){
  index = i;
  if(i > o->nConv || i<0){ //initialize to invalid values for the MVA
    pPair.SetXYZ( -999., -999.,-999.);
    pRefittedPair.SetXYZ( -999., -999.,-999.);
    vtx.SetXYZ( -999., -999.,-999.);
    eOverP = -999;
    vtxChi2 = 0;
    vtxChi2Prob = 0;
    vtxNTracks  = 0;
    vtxIsValid = 0;
    vtxMVA     = -999;
  }else{
    pPair.SetXYZ(o->pxPairConv[i],o->pyPairConv[i],o->pzPairConv[i]);
    pRefittedPair.SetXYZ(o->pxRefittedPairConv[i],o->pyRefittedPairConv[i],o->pzRefittedPairConv[i]);
    vtx.SetXYZ(o->xVtxConv[i],o->yVtxConv[i],o->zVtxConv[i]);

    p4RefittedPair.SetPtEtaPhiE(o->ptRefittedPairConv[i], 
				o->etaRefittedPairConv[i], 
				o->phiRefittedPairConv[i], 
				o->energyRefittedPairConv[i]);
    eOverP = o->eOverPRefittedPairConv[i];
    
    vtxChi2 = o->chi2VtxConv[i];
    vtxChi2Prob = o->chi2ProbVtxConv[i];
    vtxIsValid  = o->isValidVtxConv[i];
    vtxNTracks  = o->nTracksVtxConv[i];
    vtxMVA      = o->mvaOutVtxConv[i];

    trk1Dz      = o->trk1DzConv[i];
    trk1DzError = o->trk1DzErrorConv[i];  
    trk1Charge  = o->trk1ChargeConv[i];   
    trk1Algo    = o->trk1AlgoConv[i];     
    trk1D0      = o->trk1D0Conv[i];       
    trk1Pout    = o->trk1PoutConv[i];
    trk1Pin     = o->trk1PinConv[i];
    
    trk2Dz      = o->trk2DzConv[i];
    trk2DzError = o->trk2DzErrorConv[i];  
    trk2Charge  = o->trk2ChargeConv[i];   
    trk2Algo    = o->trk2AlgoConv[i];     
    trk2D0      = o->trk2D0Conv[i];       
    trk2Pout    = o->trk2PoutConv[i];
    trk2Pin     = o->trk2PinConv[i];
    
  }
}
/*
ConversionInfo VecbosConversion::getStruct(){
  ConversionInfo out = {
    index,
    pPair.X(),
    pPair.Y(),
    pPair.Z(),
    pRefittedPair.X(),
    pRefittedPair.Y(),
    pRefittedPair.Z(),
    p4RefittedPair.Pt(),
    p4RefittedPair.Eta(),
    p4RefittedPair.Phi(),
    p4RefittedPair.Energy(),
    eOverP,
    vtx.X(),
    vtx.Y(),
    vtx.Z(),
    vtxChi2,
    vtxChi2Prob,
    vtxIsValid,
    vtxNTracks,
    vtxMVA,
    trk1Dz,       
    trk1DzError,  
    trk1Charge,   
    trk1Algo,     
    trk1D0,       
    trk1Pout,
    trk1Pin,
    trk2Dz,       
    trk2DzError,  
    trk2Charge,   
    trk2Algo,     
    trk2D0,       
    trk2Pout,     
    trk2Pin
  };
  return out;
}
*/

VecbosPho::VecbosPho():
  CaloPos(0.,0.,0.)
{
  
}

VecbosPho::VecbosPho(VecbosBase* o, int i):
  CaloPos(0.,0.,0.)
{
  this->Init(o,i);
}

void VecbosPho::Init(VecbosBase* o, int i){
  if(i>=o->nPho || i<0){
    index  = -1;
    return;
  }
  SC.Init(o,o->superClusterIndexPho[i]);
  //PFSC.Init(o,SC.eta,SC.phi);
  energy = o->energyPho[i];
  eta    = o->etaPho[i];
  phi    = o->phiPho[i];
  index  = i;

  HoverE = o->hOverEPho[i];
  HTowOverE = o->hTowOverEPho[i];
  hasPixel = o->hasPixelSeedPho[i];

  dr03EcalRecHitSumEtCone = o->dr03EcalRecHitSumEtPho[i];
  dr03HcalTowerSumEtCone  = o->dr03HcalTowerSumEtPho[i];
  dr03TrkSumPtCone        = o->dr03TkSumPtPho[i];
  dr03TrkSumPtHollowCone  = o->dr03HollowTkSumPtPho[i]; 

  dr04EcalRecHitSumEtCone = o->dr04EcalRecHitSumEtPho[i];
  dr04HcalTowerSumEtCone  = o->dr04HcalTowerSumEtPho[i];
  dr04TrkSumPtCone        = o->dr04TkSumPtPho[i];
  dr04TrkSumPtHollowCone  = o->dr04HollowTkSumPtPho[i]; 

#if debugEGObject
  cout << "VecbosEGObject nPV: " << o->nPV << endl;;
#endif
  nPV=o->nPV;
  for(int iPV=0;iPV<o->nPV;iPV++){
    dr01ChargedHadronPFIso[iPV] = o->dr01ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr02ChargedHadronPFIso[iPV] = o->dr02ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr03ChargedHadronPFIso[iPV] = o->dr03ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr04ChargedHadronPFIso[iPV] = o->dr04ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr05ChargedHadronPFIso[iPV] = o->dr05ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr06ChargedHadronPFIso[iPV] = o->dr06ChargedHadronPFIsoPho[i*o->nPV+iPV];
  }
  dr01NeutralHadronPFIso  = o->dr01NeutralHadronPFIsoPho[i];
  dr02NeutralHadronPFIso  = o->dr02NeutralHadronPFIsoPho[i];
  dr03NeutralHadronPFIso  = o->dr03NeutralHadronPFIsoPho[i];
  dr04NeutralHadronPFIso  = o->dr04NeutralHadronPFIsoPho[i];
  dr05NeutralHadronPFIso  = o->dr05NeutralHadronPFIsoPho[i];
  dr06NeutralHadronPFIso  = o->dr06NeutralHadronPFIsoPho[i];

  dr01PhotonPFIso  = o->dr01PhotonPFIsoPho[i];
  dr02PhotonPFIso  = o->dr02PhotonPFIsoPho[i];
  dr03PhotonPFIso  = o->dr03PhotonPFIsoPho[i];
  dr04PhotonPFIso  = o->dr04PhotonPFIsoPho[i];
  dr05PhotonPFIso  = o->dr05PhotonPFIsoPho[i];
  dr06PhotonPFIso  = o->dr06PhotonPFIsoPho[i];

  int SCI = SC.index;
  if(SCI!=-1) CaloPos.SetXYZ(o->xPosSC[SCI],o->yPosSC[SCI],o->zPosSC[SCI]);

  this->matchConversion(o,true);
  genMatch.Init(o,-1);
};

TLorentzVector VecbosPho::p4FromVtx(TVector3 vtx,float E,bool pf){
  TVector3 dir;
  if(pf) dir = PFSC.CaloPos-vtx;
  else   dir = SC.CaloPos-vtx;
  TVector3 p   = dir.Unit()*E;
  TLorentzVector p4;
  p4.SetVectM(p,0);
  return p4;
}
#define debugConversionMatch 0
void VecbosPho::matchConversion(VecbosBase *o,bool dR){ //for some reason, the h2gglobe code doesn't always use dR
  float minDR   = 999; int iMinDR   = -1;
  float minDeta = 999; 
  float minDphi = 999; int iMinDetaDphi = -1;
  if(debugConversionMatch) cout << "nConv: " << o->nConv << endl;
  for(int i=0; i<o->nConv; i++){
    if(o->ptRefittedPairConv[i] < 1 ||
       !o->isValidVtxConv[i] ||
       o->nTracksVtxConv[i] != 2 ||
       o->chi2ProbVtxConv[i]< 1e-6) continue;  //skip bad conversions
    if(debugConversionMatch) cout << "Conversion " << i << " valid" << endl;
    float convPhi = o->phiRefittedPairConv[i];
    //convert the eta based on the z or the PV
    float convEta = etaTransformation(o->etaRefittedPairConv[i], o->zOfPVFromTracksConv[i] ); 
    if(debugConversionMatch) cout << "eta/phi: " << convEta << " / " << convPhi << endl;
    float Deta = fabs(convEta - this->SC.eta);
    float Dphi = DeltaPhi(convPhi,this->SC.phi);
    float DR   = DeltaR(convEta,this->SC.eta,convPhi,this->SC.phi);
    if(debugConversionMatch) cout << "dEta/dPhi: " << Deta << " / " << Dphi << endl;
    if(DR < minDR){
      minDR = DR;
      iMinDR = i;
    }
    if( Deta < minDeta && Dphi < minDphi){
      minDphi = Dphi; minDeta = Deta;
      iMinDetaDphi = i;
    }
  }
  int matchIndex = -1;
  if(dR){
    if(minDR < 0.1) matchIndex = iMinDR;
  }else{
    if(minDeta<0.1 && minDphi < 0.1) matchIndex = iMinDetaDphi;
  }

  conversion.Init(o,matchIndex);
}

void VecbosPho::doGenMatch(VecbosBase* o){
  const int phoID = 22;
  const float maxDR = 0.2;
  float dEoEBest = 9999;
  int indexGen = -1;
  for(int i=0;i<o->nMc;i++){  
    if(!(o->statusMc[i]==1)) continue; //require status 1 particles
    if(!(o->idMc[i] == phoID)) continue; //gen photon
    if(o->energyMc[i] < 1.) continue;
    if(DeltaR(SC.eta,o->etaMc[i],SC.phi,o->phiMc[i]) > maxDR) continue;
    float dEoE = fabs(finalEnergy-o->energyMc[i])/o->energyMc[i];
    if(dEoE > 1.) continue;

    if(dEoE < dEoEBest){
      dEoEBest = dEoE;
      indexGen = i;
    }
  }
  
  genMatch.Init(o,indexGen);
}

/*
PhoInfo VecbosPho::getStruct(){
  PhoInfo out={
    index,
    energy,
    eta,
    phi,
    correctedEnergy,
    correctedEnergyError,
    scaledEnergy,
    smearedEnergy,
    dEoE,
    HoverE,
    hasPixel,
    CaloPos.X(),
    CaloPos.Y(),
    CaloPos.Z(),
    dr03EcalRecHitSumEtCone,
    dr03HcalTowerSumEtCone,
    dr03TrkSumPtCone,
    dr03TrkSumPtHollowCone,
    dr04EcalRecHitSumEtCone,
    dr04HcalTowerSumEtCone,
    dr04TrkSumPtCone,
    dr04TrkSumPtHollowCone,
    conversion->index,
    SC.index
  };
  return out;
}
*/

VecbosEle::VecbosEle(){}

VecbosEle::VecbosEle(VecbosBase* o,int i){
  this->Init(o,i);
}
void VecbosEle::Init(VecbosBase* o,int i){
  if(i>o->nEle || i < 0){
    index = -1;
    return;
  }
  SC.Init(o,o->superClusterIndexPho[i]);
  if(SC.index==-1){  // an electron without a supercluster??
    index = -1; 
    return;
  }
  energy = o->energyEle[i];
  correctedEnergy = 0;

  eta    = o->etaEle[i];
  phi    = o->phiEle[i];
  pt     = TMath::Sqrt(TMath::Power(o->pxEle[i],2)+TMath::Power(o->pyEle[i],2));

  charge = o->chargeEle[i];

  esEnergy = SC.esEnergy;
  HoverE   = o->hOverEPho[i];

  vtxX = o->vertexXEle[i];
  vtxY = o->vertexYEle[i];
  vtxZ = o->vertexZEle[i];

  EOverP = o->eSuperClusterOverPEle[i];

  dr03ChargedHadronPFIso = o->pfCandChargedIso03Ele[i];
  dr03NeutralHadronPFIso = o->pfCandNeutralIso03Ele[i];
  dr03PhotonPFIso        = o->pfCandPhotonIso03Ele[i];

  dr04ChargedHadronPFIso = o->pfCandChargedIso04Ele[i];
  dr04NeutralHadronPFIso = o->pfCandNeutralIso04Ele[i];
  dr04PhotonPFIso        = o->pfCandPhotonIso04Ele[i];

  dr03TkSumPt          = o->dr03TkSumPtEle[i];
  dr03EcalRecHitSumEt  = o->dr03EcalRecHitSumEtEle[i];
  dr03HcalTowerSumEt   = o->dr03HcalTowerSumEtEle[i];

  dr04TkSumPt          = o->dr04TkSumPtEle[i];
  dr04EcalRecHitSumEt  = o->dr04EcalRecHitSumEtEle[i];
  dr04HcalTowerSumEt   = o->dr04HcalTowerSumEtEle[i];

  dEtaSCTrackAtCalo = o->deltaEtaAtCaloEle[i];
  dPhiSCTrackAtCalo = o->deltaPhiAtCaloEle[i];

  dEtaSCTrackAtVtx = o->deltaEtaAtVtxEle[i];
  dPhiSCTrackAtVtx = o->deltaPhiAtVtxEle[i];

  idMVA  = o->mvaidnontrigEle[i];

  int iGsfTrack = o->gsfTrackIndexEle[i];
  if(iGsfTrack < 0 || iGsfTrack >= o->nGsfTrack){
    d0Track=999;
    dzTrack=999;
    expInnerLayersHits=-999;
  }else{
    d0Track= o->d0GsfTrack[iGsfTrack];
    dzTrack= o->dzGsfTrack[iGsfTrack];
    expInnerLayersHits = o->expInnerLayersGsfTrack[iGsfTrack];
  }

  int tmp = o->recoFlagsEle[i];
  //extract the flag information into booleans
  isTrackerDriven = tmp & 1; 
  isEcalDriven    = (tmp >> 1) & 1;

  hasMatchedConversion = o->hasMatchedConversionEle[i];

  genMatch.Init(o,-1);
};

void VecbosEle::doGenMatch(VecbosBase* o){
  const int eleID = 11;
  const float maxDR = 0.2;
  float dEoEBest = 9999;
  int indexGen = -1;
  for(int i=0;i<o->nMc;i++){  
    if(!o->statusMc[i]==1) continue; //require status 1 particles
    if(!abs(o->idMc[i]) == eleID) continue; //gen photon
    if(o->energyMc[i] < 1.) continue;
    if(DeltaR(SC.eta,o->etaMc[i],SC.phi,o->phiMc[i]) > maxDR) continue;
    float dEoE = fabs(energy-o->energyMc[i])/o->energyMc[i];
    if(dEoE > 1.) continue;
    if(dEoE < dEoEBest){
      dEoEBest = dEoE;
      indexGen = i;
    }
  }
  genMatch.Init(o,indexGen);
}

/*
EleInfo VecbosEle::getStruct(){
  EleInfo out={
  index,
  energy,
  eta,
  phi,
  esEnergy,
  HoverE,
  isEcalDriven,
  isTrackerDriven,
  SC.index,
  CaloPos.X(),  
  CaloPos.Y(),  
  CaloPos.Z()  
  };
  return out;
}
*/

VecbosMu::VecbosMu(){}

VecbosMu::VecbosMu(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosMu::Init(VecbosBase* o, int i){
  if(i<0 || i > o->nMuon){
    index = -1;
  }else{
    index = i;            
    energy = o->energyMuon[i];         
    pt = TMath::Sqrt( TMath::Power(o->pxMuon[i],2) + TMath::Power(o->pyMuon[i],2) );
    eta = o->etaMuon[i];
    phi = o->phiMuon[i];
    //p4.SetPtEtaPhiE(pt,eta,phi,energy);
    charge = o->chargeMuon[i];
    combinedIso = (o->emEt03Muon[i] + o->hadEt03Muon[i] + o->sumPt03Muon[i] - o->rhoFastjet * TMath::Pi()*0.3*0.3)/pt;;
    emIso = o->emEt03Muon[i];
    hadIso = o->hadEt03Muon[i];
    trkIso = o->sumPt03Muon[i];
    Utils AnalysisUtilities;
    isGlobalMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::AllGlobalMuons);
    isTrackerMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::AllTrackerMuons);
    isPromptMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::GlobalMuonPromptTight);
    int iTrack = o->trackIndexMuon[i];
    if(iTrack<0 || iTrack > o->nTrack) {
      nTrackHits=0;
      nPixelHits=0;
      trackImpactPar=0;
    }else{
      nTrackHits = o->numberOfValidStripTIBHitsTrack[iTrack]
	+ o->numberOfValidStripTIDHitsTrack[iTrack]
	+ o->numberOfValidStripTOBHitsTrack[iTrack]
	+ o->numberOfValidStripTECHitsTrack[iTrack];
      nPixelHits = o->numberOfValidPixelBarrelHitsTrack[iTrack] + o->numberOfValidPixelEndcapHitsTrack[iTrack];
      trackImpactPar = fabs(o->transvImpactParTrack[iTrack]);
    }
    isLooseMuon = true;
    if(!isGlobalMuon || nTrackHits<=10) isLooseMuon=false;
    isTightMuon = isLooseMuon;
    if(!isTrackerMuon || !isPromptMuon || combinedIso >=0.15 || nPixelHits == 0 || trackImpactPar >=0.2) isTightMuon = false;
  }

  genMatch.Init(o,-1);
}

void VecbosMu::doGenMatch(VecbosBase* o){
  const int muID = 13;
  const float maxDR = 0.2;
  float dEoEBest = 9999;
  int indexGen = -1;
  for(int i=0;i<o->nMc;i++){  
    if(!o->statusMc[i]==1) continue; //require status 1 particles
    if(!abs(o->idMc[i]) == muID) continue; //gen mu
    if(o->energyMc[i] < 1.) continue;
    if(DeltaR(eta,o->etaMc[i],phi,o->phiMc[i]) > maxDR) continue;
    float dEoE = fabs(energy-o->energyMc[i])/o->energyMc[i];
    if(dEoE > 0.5) continue;
    if(dEoE < dEoEBest){
      dEoEBest = dEoE;
      indexGen = i;
    }
  }
  genMatch.Init(o,indexGen);
}



VecbosGen::VecbosGen(){}

VecbosGen::VecbosGen(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosGen::Init(VecbosBase* o, int i){
  if(i < 0 || i >= o->nMc){
    index = -1;
    return;
  }

  index = i;
  pt = o->pMc[i]/cosh(o->etaMc[i]);
  energy = o->energyMc[i];
  eta = o->etaMc[i];
  phi = o->phiMc[i];
  mass = TMath::Sqrt( TMath::Power(o->energyMc[i],2)-TMath::Power(o->pMc[i],2) );
  
  //Vtx.SetXYZ(o->vxMc[i],o->vyMc[i],o->vzMc[i]);
  Vx = o->vxMc[i];
  Vy = o->vyMc[i];
  Vz = o->vzMc[i];

  status = o->statusMc[i];
  id     = o->idMc[i];
  indexMother = o->mothMc[i];
  
  if(indexMother >=0 && indexMother < o->nMc){
    statusMother = o->statusMc[indexMother];
    idMother     = o->idMc[indexMother];
    
    if( id == idMother ){ //mother is the same particle: probably documentation line.
      int indexMo_tmp = indexMother;
      
      for(int i=0;i<indexMother;i++) // find the true mother
  	{
	  indexMother = o->mothMc[indexMo_tmp];
  	  idMother = o->idMc[indexMother];
	  indexMo_tmp = indexMother;

  	  if ( idMother != id ) break;
  	}
    }
  }else{
    statusMother = -1;
    idMother     =  0;
  }
}

/*
TLorentzVector VecbosGen::getP4(){
  TLorentzVector p4;
  p4.SetPtEtaPhiE(pt,eta,phi,energy);
  return p4;
}
*/


//jets 

VecbosJet::VecbosJet(){}

VecbosJet::VecbosJet(VecbosBase* o, int i,VecbosJet::JetType jtype = PFPUcorr){
  this->Init(o,i,jtype);
}

void VecbosJet::Init(VecbosBase* o, int i,VecbosJet::JetType jtype = PFPUcorr){
  type = jtype;
  switch(type){
  case PFPUcorr:
    if(i < 0 || i > o->nAK5PFPUcorrJet){
      index = -1;
      return;
    }
    index = i;
   energy = o->energyAK5PFPUcorrJet[i];
   uncorrEnergy = o->uncorrEnergyAK5PFPUcorrJet[i];
   uncorrpx = o->uncorrpxAK5PFPUcorrJet[i];
   uncorrpy = o->uncorrpyAK5PFPUcorrJet[i];
   uncorrpz = o->uncorrpzAK5PFPUcorrJet[i];
   eta = o->etaAK5PFPUcorrJet[i];
   phi = o->phiAK5PFPUcorrJet[i];
   pt = TMath::Sqrt(TMath::Power(o->pxAK5PFPUcorrJet[i],2) + TMath::Power(o->pyAK5PFPUcorrJet[i],2));
   
   charge = o->chargeAK5PFPUcorrJet[i];

   vtxX = o->vertexXAK5PFPUcorrJet[i];
   vtxY = o->vertexYAK5PFPUcorrJet[i];
   vtxZ = o->vertexZAK5PFPUcorrJet[i];

   area =  o->areaAK5PFPUcorrJet[i];
   chargedHadronFraction = o->chargedHadronEnergyAK5PFPUcorrJet[i];
   neutralHadronFraction = o->neutralHadronEnergyAK5PFPUcorrJet[i];

   jetIdMva = o->jetIdMvaAK5PFPUcorrJet[i];
  
   betaStar = o->betastarAK5PFPUcorrJet[i];
   betaStarIdMVA = o->betaIdMvaAK5PFPUcorrJet[i];
   betaStarClassicIdMVA = o->betastarclassicIdMvaAK5PFPUcorrJet[i];
   rmsCands = o->rmsCandAK5PFPUcorrJet[i];
   rmsCandsHand = o->rmsCandsHandAK5PFPUcorrJet[i];

   combinedSecondaryVertex = o->combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[i];
   simpleSecondaryVertexHighPur = o->simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet[i];
   simpleSecondaryVertexHighEff = o->simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[i];
   break;

  case PFNoPU:
        if(i < 0 || i > o->nAK5PFNoPUJet){
      index = -1;
      return;
    }
    index = i;
   energy = o->energyAK5PFNoPUJet[i];
   uncorrEnergy = o->uncorrEnergyAK5PFNoPUJet[i];
   uncorrpx = o->uncorrpxAK5PFNoPUJet[i];
   uncorrpy = o->uncorrpyAK5PFNoPUJet[i];
   uncorrpz = o->uncorrpzAK5PFNoPUJet[i];
   eta = o->etaAK5PFNoPUJet[i];
   phi = o->phiAK5PFNoPUJet[i];
   pt = TMath::Sqrt(TMath::Power(o->pxAK5PFNoPUJet[i],2) + TMath::Power(o->pyAK5PFNoPUJet[i],2));
   
   charge = o->chargeAK5PFNoPUJet[i];

   vtxX = o->vertexXAK5PFNoPUJet[i];
   vtxY = o->vertexYAK5PFNoPUJet[i];
   vtxZ = o->vertexZAK5PFNoPUJet[i];

   area =  o->areaAK5PFNoPUJet[i];
   chargedHadronFraction = o->chargedHadronEnergyAK5PFNoPUJet[i];
   neutralHadronFraction = o->neutralHadronEnergyAK5PFNoPUJet[i];

   jetIdMva = o->jetIdMvaAK5PFNoPUJet[i];
  
   betaStar = o->betastarAK5PFNoPUJet[i];
   betaStarIdMVA = o->betaIdMvaAK5PFNoPUJet[i];
   betaStarClassicIdMVA = o->betastarclassicIdMvaAK5PFNoPUJet[i];
   rmsCands = o->rmsCandAK5PFNoPUJet[i];
   rmsCandsHand = o->rmsCandsHandAK5PFNoPUJet[i];

   combinedSecondaryVertex = o->combinedSecondaryVertexBJetTagsAK5PFNoPUJet[i];
   simpleSecondaryVertexHighPur = o->simpleSecondaryVertexHighPurBJetTagsAK5PFNoPUJet[i];
   simpleSecondaryVertexHighEff = o->simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet[i];
   break;
 
  default:
   index = -1;
   break;
  }

}


TLorentzVector VecbosJet::getP4(){
  TLorentzVector p4;
  p4.SetPtEtaPhiE(pt,eta,phi,energy);
  return p4;
}
