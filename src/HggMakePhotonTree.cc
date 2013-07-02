#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

#include "HggPhysUtils.cc"
#include "ecalGap.cc"
#include <iostream>

#include "HggMakePhotonTree.hh"

#define debugPhotonTree 0
using namespace std;
HggMakePhotonTree::HggMakePhotonTree(TTree* t, string outputFile,bool data):
  VecbosBase(t),
  usePF(false)
{

  outFName = outputFile;
  isRealData = data;
}

void HggMakePhotonTree::Loop(int start=0,int stop=-1){
  this->init();
  Long64_t jentry = start-1;
  while(fChain->GetEntry(++jentry)){
    if(jentry == stop) break;
    if(jentry %1000==0)cout << "Processing Entry: " << jentry << endl;
    //do this for each photon
    //    cout << nPho << endl;
    for(int iPho=0;iPho<nPho;iPho++){      
      PhoIndex = iPho;
      SCPhoIndex = superClusterIndexPho[iPho];
      if(SCPhoIndex <0 || SCPhoIndex > nSC) continue; //no super cluster for a photon??? thats weird ...
      if(debugPhotonTree) cout << PhoIndex << "  " << SCPhoIndex << endl;      
      GenPhoIndex = genMatchPho();
      if(debugPhotonTree) cout << GenPhoIndex << endl;	    
      truephtmatched = (GenPhoIndex != -1); 
      phoPt = TMath::Sqrt(pxPho[iPho]*pxPho[iPho] + pyPho[iPho]*pyPho[iPho]);
      if(phoPt < 28.) continue;
      ConvPhoIndex = this->convMatchPho();
      if(debugPhotonTree) cout << ConvPhoIndex << endl;
      if(debugPhotonTree) cout << "Gen Info" << endl;
      this->fillGenInfo();      
      if(debugPhotonTree) cout << "Gen Info" << endl;
      this->fillGapVariables();
      if(debugPhotonTree) cout << "Gap Variables" << endl;

      this->fillEnergyVariables();
      if(debugPhotonTree) cout << "Fill Energy Variables" << endl;
      this->fillIsoVariables();
      if(debugPhotonTree) cout << "Fill IsoVariables" << endl;
      this->fillConvVariables();
      if(debugPhotonTree) cout << "Fill Conv Variables" << endl;
      nVertex = nPV;
      rho  = rhoFastjet;
      haspromptele=0;
      evtNumber = eventNumber;
      if(debugPhotonTree) cout << "Looping Over Electrons" << endl;
      for(int iEle=0;iEle<nEle;iEle++){
	//copied from  ConversionTools::hasMatchedPromptElectron
	if(superClusterIndexEle[iEle]!=SCPhoIndex) continue;
	if(hasMatchedConversionEle[iEle]) continue;
	if(gsfTrackIndexEle[iEle]<0 || gsfTrackIndexEle[iEle] >=nGsfTrack) continue;
	if(expInnerLayersGsfTrack[gsfTrackIndexEle[iEle]]>0) continue;
	haspromptele=1;
	break;
      }

      tree->Fill();
   }    
  }//end main loop
  if(debugPhotonTree) cout << "DONE: Writing" << endl;
  TFile *f = new TFile(outFName.c_str(),"RECREATE");
  tree->Write();
  f->Close();

}

void HggMakePhotonTree::init(){
  isRealData=false;
  ecalGeometry = loadecalGapCoordinates(isRealData);
  this->setupBranches();
}

void HggMakePhotonTree::fillGenInfo(){
  if(GenPhoIndex==-1){
    etrue = -99;
    etatrue = -99;
    phitrue = -99;
    pidphtmom = -99;
    etaecaltrue = -99;
    phiecaltrue = -99;
  }
  etrue = energyMc[GenPhoIndex];
  etatrue = etaMc[GenPhoIndex];
  phitrue  = phiMc[GenPhoIndex];
  if(mothMc[GenPhoIndex] >=0)
    pidphtmom = idMc[mothMc[GenPhoIndex]];
  else
    pidphtmom = 0;
  etaecaltrue =  ecalEta(etaMc[GenPhoIndex],vzMc[GenPhoIndex],
			 sqrt(vxMc[GenPhoIndex]*vxMc[GenPhoIndex]+vyMc[GenPhoIndex]*vyMc[GenPhoIndex]));
  phiecaltrue =  ecalPhi(phiMc[GenPhoIndex],vxMc[GenPhoIndex],vyMc[GenPhoIndex]);
}

void HggMakePhotonTree::fillGapVariables(){
  THIS_ECAL_GEO geo = getGapCoordinates(ecalGeometry,etaSC[SCPhoIndex],phiSC[SCPhoIndex]);
  etaCGap = geo._aC;
  phiCGap = geo._bC;
  etaSGap = geo._aS;
  phiSGap = geo._bS;
  etaMGap = geo._aM;
  phiMGap = geo._bM;

}

void HggMakePhotonTree::fillIsoVariables(){
  ecalisodr03 = dr03EcalRecHitSumEtPho[PhoIndex] - rho*rhoFac;         
  ecalisodr04 = dr04EcalRecHitSumEtPho[PhoIndex] - rho*rhoFac;         
  hcalisodr03 = dr03HcalTowerSumEtPho[PhoIndex]  - rho*rhoFac;           
  hcalisodr04 = dr04HcalTowerSumEtPho[PhoIndex]  - rho*rhoFac;           
  hoe = hOverEPho[PhoIndex];

  float worstIsoDR04=0; 
  int   indexWorstIsoDR04=-1;
  int   indexGenMatchVtx=-1;
  TVector3 genVtx(vxMc[GenPhoIndex],vyMc[GenPhoIndex],vzMc[GenPhoIndex]);
  float minVtxD=999;
  for(int iVtx=0;iVtx<nPV;iVtx++){
    float thisIsoDR03 = computeTrackIso(PhoIndex,iVtx,0.,0.3,0.2,0.0,1.0,0.1);
    float thisIsoDR04 = computeTrackIso(PhoIndex,iVtx,0.,0.4,0.2,0.0,1.0,0.1);
    if(thisIsoDR04 > worstIsoDR04){
      worstIsoDR04 = thisIsoDR04;
      indexWorstIsoDR04 = iVtx;
    }
    TVector3 thisVtx(PVxPV[iVtx],PVyPV[iVtx],PVzPV[iVtx]);
    if( (genVtx-thisVtx).Mag() < minVtxD){
      minVtxD = (genVtx-thisVtx).Mag();
      indexGenMatchVtx = iVtx;
    }
  }
  trkisodr03vtxbestGenMatched = computeTrackIso(PhoIndex,indexGenMatchVtx,0.,0.3,0.2,0.0,1.0,0.1);
  isosumoet=(trkisodr03vtxbestGenMatched +ecalisodr03  + hcalisodr04 -rho*rhoFac)*50./phoPt;
  trkisooet=(trkisodr03vtxbestGenMatched)*50./phoPt;
  TLorentzVector phop4_bad = p4FromVtx(superClusterIndexPho[PhoIndex],indexWorstIsoDR04,energyPho[PhoIndex]);
  isosumoetbad=(trkisodr04vtxWorst + ecalisodr04  + hcalisodr04 -rho*rhoFacBad)*50./phop4_bad.Et();
}

void HggMakePhotonTree::fillConvVariables(){
  if(ConvPhoIndex==-1){
    convp = -1;              
    convpt = -1;             
    convpttrk1 = -1;         
    convpttrk2 = -1;         
    deltaphi_convsc = -1;    
    deltaeta_convsc = -1;    
    convz = -1;              
    convrho = -1;  
  }else{ // we have a conversions!!
    convp  = energyRefittedPairConv[ConvPhoIndex];
    convpt = ptRefittedPairConv[ConvPhoIndex];
    float eta = etaTransformation(etaRefittedPairConv[ConvPhoIndex],zOfPVFromTracksConv[ConvPhoIndex]);
    float phi = ecalPhi(phiRefittedPairConv[ConvPhoIndex],xVtxConv[ConvPhoIndex],yVtxConv[ConvPhoIndex]);
    convpttrk1 = trk1PtConv[ConvPhoIndex];
    convpttrk2 = trk2PtConv[ConvPhoIndex];
    deltaphi_convsc = DeltaPhi(phiSC[SCPhoIndex],phi);
    deltaeta_convsc = fabs(etaSC[SCPhoIndex]-eta);
    convz = zVtxConv[ConvPhoIndex];
    convrho = TMath::Sqrt( TMath::Power(xVtxConv[ConvPhoIndex],2) + TMath::Power(yVtxConv[ConvPhoIndex],2) );
  }
}

struct sort_pred {
  bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right) {
    return left.second < right.second; 
  }
}; 
void HggMakePhotonTree::fillEnergyVariables(){
  escraw = rawEnergySC[SCPhoIndex];
  e5x5   = e5x5SC[SCPhoIndex];
  e2x2   = e2x2SC[SCPhoIndex];
  e3x3   = e3x3SC[SCPhoIndex];
  e3x1   = e3x1SC[SCPhoIndex];
  e1x3   = e1x3SC[SCPhoIndex];
  e4x4   = e4x4SC[SCPhoIndex];
  e2x5   = e2x5MaxSC[SCPhoIndex];
  e2x5right  = e2x5RightSC[SCPhoIndex];
  e2x5left   = e2x5LeftSC[SCPhoIndex];
  e2x5top    = e2x5TopSC[SCPhoIndex];
  e2x5bottom = e2x5BottomSC[SCPhoIndex];
  e2x5max    = e2x5MaxSC[SCPhoIndex];
  eps    = esEnergySC[SCPhoIndex];
  esc    = energySC[SCPhoIndex];
  etasc  = etaSC[SCPhoIndex];
  phisc  = phiSC[SCPhoIndex];
  scetawidth = etaWidthSC[SCPhoIndex];
  scphiwidth = phiWidthSC[SCPhoIndex];
  sigietaieta    = sqrt(covIEtaIEtaSC[SCPhoIndex]);
  sigcovietaiphi = covIEtaIPhiSC[SCPhoIndex];
  sigiphiiphi    = sqrt(covIPhiIPhiSC[SCPhoIndex]);
  emax   = eMaxSC[SCPhoIndex];
  scnbc  = nBCSC[SCPhoIndex];
  scncrystal = nCrystalsSC[SCPhoIndex];
  etop     = eTopSC[SCPhoIndex];
  ebottom  = eBottomSC[SCPhoIndex];
  eleft    = eLeftSC[SCPhoIndex];
  eright   = eRightSC[SCPhoIndex];
  e2max    = e2ndSC[SCPhoIndex];

  scbcfmax = 0;
  int seedIndex = -1;
  std::vector<std::pair<int,float> > bcs;
  for(int iBC=0;iBC<nBC;iBC++){
    if(indexSCBC[iBC]!=SCPhoIndex) continue;
    bcs.push_back(std::pair<int,float>(iBC,energyBC[iBC]));
    if(energyBC[iBC] > scbcfmax){
      scbcfmax = energyBC[iBC];
      seedIndex = iBC;
    }
  }
  std::sort(bcs.begin(),bcs.end(),sort_pred());

  scbcfmax/=rawEnergySC[SCPhoIndex];

  epht   = energyPho[PhoIndex];
  etapht = etaPho[PhoIndex];
  phipht = phiPho[PhoIndex];
  r9 = e3x3SC[SCPhoIndex]/rawEnergySC[SCPhoIndex];


  seedieta = iEtaBC[seedIndex];
  seediphi = iPhiBC[seedIndex];
  bcseedetacry = etaCrystalBC[seedIndex];
  bcseedphicry = phiCrystalBC[seedIndex];

  //basic clusters
  int i1 = bcs.at(0).first;
  bce2nd = e2ndBC[i1];
  dbceta = etaSC[SCPhoIndex] - etaBC[i1];
  dbcphi = DeltaPhi(phiSC[SCPhoIndex],phiBC[i1]);
  bce    = energyBC[i1];

  bool b2 = (bcs.size() > 1);
  int i2  = (b2 ? bcs.at(1).first:0);

  bc2e     = (b2 ? energyBC[i2] : -99);
  bc2emax  = (b2 ? eMaxBC[i2] : -99);
  bc2e2nd  = (b2 ? e2ndBC[i2] : -99);
  bc2etop     = (b2 ? eTopBC[i2] : -99);
  bc2ebottom  = (b2 ? eBottomBC[i2] : -99);
  bc2eleft    = (b2 ? eLeftBC[i2] : -99);
  bc2eright   = (b2 ? eRightBC[i2] : -99);
  bc2sigietaieta = (b2 ? sqrt(covIEtaIEtaBC[i2]) : -99);
  bc2sigiphiiphi = (b2 ? sqrt(covIPhiIPhiBC[i2]) : -99);
  bc2sigietaiphi = (b2 ? covIEtaIPhiBC[i2] : -99);
  bc2e3x3   = (b2 ? e3x3BC[i2] : -99);
  bc2e5x5   = (b2 ? e5x5BC[i2] : -99);
  dbc2eta   =  (b2 ? etaSC[SCPhoIndex] - etaBC[i2] : -99);
  dbc2phi   =  (b2 ? DeltaPhi(phiSC[SCPhoIndex],phiBC[i2]) : -99);

  bc2ieta   = (b2 ? iEtaBC[i2] : -99);
  bc2iphi   = (b2 ? iPhiBC[i2] : -99);
  bc2etacry = (b2 ? etaCrystalBC[i2] : -99);
  bc2phicry = (b2 ? phiCrystalBC[i2] : -99);

  dbc2bceta = (b2 ? etaBC[i1]-etaBC[i2] : -99);
  dbc2bcphi = (b2 ? DeltaPhi(phiBC[i1],phiBC[i2]) : -99);
  
  bool blast = (bcs.size() > 2);
  int  ilast = (blast ? bcs.back().first : 0);
  bclaste     = (blast ? energyBC[ilast] : -99);
  bclaste3x3   = (blast ? e3x3BC[ilast] : -99);
  bclaste5x5   = (blast ? e5x5BC[ilast] : -99);
  bclastsigietaieta = (blast ? sqrt(covIEtaIEtaBC[ilast]) : -99);
  bclastsigiphiiphi = (blast ? sqrt(covIPhiIPhiBC[ilast]) : -99);
  bclastsigietaiphi = (blast ? covIEtaIPhiBC[ilast] : -99);
  dbclasteta   =  (blast ? etaSC[SCPhoIndex] - etaBC[ilast] : -99);
  dbclastphi   =  (blast ? DeltaPhi(phiSC[SCPhoIndex],phiBC[ilast]) : -99);
  
  bool b2last = (bcs.size() > 3);
  int  i2last = (b2last ? bcs.at(bcs.size()-2).first : 0);
  bclast2e     = (b2last ? energyBC[i2last] : -99);
  bclast2e3x3   = (b2last ? e3x3BC[i2last] : -99);
  bclast2e5x5   = (b2last ? e5x5BC[i2last] : -99);
  bclast2sigietaieta = (b2last ? sqrt(covIEtaIEtaBC[i2last]) : -99);
  bclast2sigiphiiphi = (b2last ? sqrt(covIPhiIPhiBC[i2last]) : -99);
  bclast2sigietaiphi = (b2last ? covIEtaIPhiBC[i2last] : -99);
  dbclast2eta   =  (b2last ? etaSC[SCPhoIndex] - etaBC[i2last] : -99);
  dbclast2phi   =  (b2last ? DeltaPhi(phiSC[SCPhoIndex],phiBC[i2last]) : -99);
  
}

int HggMakePhotonTree::genMatchPho(){
  float maxDR=0.08;
  if(fabs(etaSC[SCPhoIndex]) < 1.48) maxDR = 0.05;
  float minDR=999; int minInd=-1;
  for(int i=0;i<nMc;i++){
    if(statusMc[i]!=1) continue;
    if(idMc[i]!=22) continue; //require photon
    float thisDR = DeltaR(etaMc[i],etaPho[PhoIndex],phiMc[i],phiPho[PhoIndex]);
    if(thisDR<minDR && thisDR<maxDR){
      minDR = thisDR; minInd = i;
    }
  }
  return minInd;
}

int HggMakePhotonTree::convMatchPho(){
  int index = SCPhoIndex;
  const float maxdEta = 0.1;
  const float maxdPhi = 0.1;
  float mindEta = 9999; float mindPhi = 9999;
  int minIndex = -1;
  for(int iConv=0;iConv<nConv;iConv++){
    if(ptRefittedPairConv[iConv] < 1.) continue;

    if(isValidVtxConv[iConv] == 0
       || nTracksVtxConv[iConv] !=2
       || chi2ProbVtxConv[iConv] <0.0005) continue;

    float convETA = etaTransformation(etaRefittedPairConv[iConv],zOfPVFromTracksConv[iConv]);
    float dEta = TMath::Power(etaSC[index]-convETA,2);
    float dPhi = TMath::Power(DeltaPhi(phiSC[index],phiRefittedPairConv[iConv]),2);

    if(dEta < mindEta && dPhi < mindPhi && dEta < maxdEta && dPhi < maxdPhi){
      mindEta = dEta; mindPhi = dPhi;
      minIndex = iConv;
    }
  }
  return minIndex;
}

TLorentzVector HggMakePhotonTree::p4FromVtx(int sc,int vtx,float E){
  TVector3 caloPos(xPosSC[sc],yPosSC[sc],zPosSC[sc]);
  TVector3 vtxPos(PVxPV[vtx],PVyPV[vtx],PVzPV[vtx]);
  TVector3 dir = caloPos-vtx;
  TVector3 p = dir.Unit()*E;
  TLorentzVector p4(p.x(),p.y(),p.z(),E);
  return p4;

}

float HggMakePhotonTree::computeTrackIso(int iPho, int iVtx,
					 float ptMin,
					 float outerCone,
					 float innerCone,
					 float etaStripHalfWidth,
					 float dzMax,
					 float dxyMax){

  float vX=0,vY=0,vZ=0;
  if(iVtx >= -0 && iVtx < nPV){
    vX = PVxPV[iVtx];
    vY = PVyPV[iVtx];
    vZ = PVzPV[iVtx];
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
    double dEta = fabs(etaPho[iPho] - trackP.Eta());
    double dR   = DeltaR(etaPho[iPho],trackP.Eta(),phiPho[iPho],trackP.Phi());
    if( dR < outerCone && dR >= innerCone && dEta >= etaStripHalfWidth) SumTrackPt+=ptTrack;
  }
  return SumTrackPt;
}
void HggMakePhotonTree::setupBranches(){
  tree = new TTree("Analysis","selected photon pair");

  tree->Branch("etaCGap",&etaCGap,"etaCGap/F");
  tree->Branch("phiCGap",&phiCGap,"phiCGap/F");
  tree->Branch("etaSGap",&etaSGap,"etaSGap/F");
  tree->Branch("phiSGap",&phiSGap,"phiSGap/F");
  tree->Branch("etaMGap",&etaMGap,"etaMGap/F");
  tree->Branch("phiMGap",&phiMGap,"phiMGap/F");
  
  tree->Branch("isosumoet",&isosumoet,"isosumoet/F");
  tree->Branch("isosumoetbad",&isosumoetbad,"isosumoetbad/F");
  tree->Branch("trkisooet",&trkisooet,"trkisooet/F");
  tree->Branch("hoe",&hoe,"hoe/F");
  tree->Branch("r9",&r9,"r9/F");
  tree->Branch("haspromptele",&haspromptele,"haspromptele/I");
  
  tree->Branch("ecalisodr03",&ecalisodr03,"ecalisodr03/F");
  tree->Branch("ecalisodr04",&ecalisodr04,"ecalisodr04/F");
  tree->Branch("hcalisodr03",&hcalisodr03,"hcalisodr03/F");
  tree->Branch("hcalisodr04",&hcalisodr04,"hcalisodr04/F");
  tree->Branch("trkisodr03vtxbestGenMatched",&trkisodr03vtxbestGenMatched,"trkisodr03vtxbestGenMatched/F");
  tree->Branch("trkisodr04vtxWorst",&trkisodr04vtxWorst,"trkisodr04vtxWorst/F");
  
  
  tree->Branch("etrue",&etrue,"etrue/F");
  tree->Branch("etatrue",&etatrue,"etatrue/F");
  tree->Branch("phitrue",&phitrue,"phitrue/F");
  tree->Branch("etaecaltrue",&etaecaltrue,"etaecaltrue/F");
  tree->Branch("phiecaltrue",&phiecaltrue,"phiecaltrue/F");
  
  
  tree->Branch("escraw",&escraw,"escraw/F");
  tree->Branch("e5x5",&e5x5,"e5x5/F");
  tree->Branch("eps",&eps,"eps/F");
  tree->Branch("esc",&esc,"esc/F");
  tree->Branch("epht",&epht,"epht/F");
  tree->Branch("etasc",&etasc,"etasc/F");
  tree->Branch("phisc",&phisc,"phisc/F");

  tree->Branch("e2x2",&e2x2,"e2x2/F");
  tree->Branch("e3x3",&e3x3,"e3x3/F");
  tree->Branch("e1x3",&e1x3,"e1x3/F");
  tree->Branch("e3x1",&e3x1,"e3x1/F");
  tree->Branch("e4x4",&e4x4,"e4x4/F");
  tree->Branch("e2x5",&e2x5,"e2x5/F");
  tree->Branch("e2x5right",&e2x5right,"e2x5right/F");
  tree->Branch("e2x5left",&e2x5left,"e2x5left/F");
  tree->Branch("e2x5top",&e2x5top,"e2x5top/F");
  tree->Branch("e2x5bottom",&e2x5bottom,"e2x5bottom/F");
  tree->Branch("e2x5max",&e2x5max,"e2x5max/F");
  

  tree->Branch("etapht",&etapht,"etapht/F");
  tree->Branch("phipht",&phipht,"phipht/F");
  
  tree->Branch("scetawidth",&scetawidth,"scetawidth/F");
  tree->Branch("scphiwidth",&scphiwidth,"scphiwidth/F");
  tree->Branch("sigietaieta",&sigietaieta,"sigietaieta/F");
  tree->Branch("sigiphiiphi",&sigiphiiphi,"sigiphiiphi/F");
  tree->Branch("sigcovietaiphi",&sigcovietaiphi,"sigcovietaiphi/F");
  tree->Branch("emax",&emax,"emax/F");
  tree->Branch("scnbc",&scnbc,"scnbc/I");
  tree->Branch("scncrystal",&scncrystal,"scncrystal/I");
  tree->Branch("scbcfmax",&scbcfmax,"scbcfmax/F");
  tree->Branch("seedieta",&seedieta,"seedieta/I");
  tree->Branch("seediphi",&seediphi,"seediphi/I");
  tree->Branch("eleft",&eleft,"eleft/F");
  tree->Branch("e2max",&e2max,"e2max/F");
  tree->Branch("eright",&eright,"eright/F");
  tree->Branch("etop",&etop,"etop/F");
  tree->Branch("ebottom",&ebottom,"ebottom/F");
  
  tree->Branch("convp",&convp,"convp/F");
  tree->Branch("convpttrk1",&convpttrk1,"convpttrk1/F");
  tree->Branch("convpttrk2",&convpttrk2,"convpttrk2/F");
  tree->Branch("convpt",&convpt,"convpt/F");
  tree->Branch("deltaphi_convsc",&deltaphi_convsc,"deltaphi_convsc/F");
  tree->Branch("deltaeta_convsc",&deltaeta_convsc,"deltaeta_convsc/F");
  tree->Branch("truephtmatched",&truephtmatched,"truephtmatched/I");
  tree->Branch("pidphtmom",&pidphtmom,"pidphtmom/I");
  
  tree->Branch("convrho",&convrho,"convrho/F");
  tree->Branch("convz",&convz,"convz/F");
  
  
  tree->Branch("rho",&rho,"rho/F");
  tree->Branch("nVertex",&nVertex,"nVertex/I");
  
  tree->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
  tree->Branch("runNumber",&runNumber,"runNumber/I");
  tree->Branch("evtNumber",&evtNumber,"evtNumber/I");
  tree->Branch("isRealData",&isRealData,"isRealData/I");
  
  tree->Branch("dbc2bceta",&dbc2bceta,"dbc2bceta/F");
  tree->Branch("dbc2bcphi",&dbc2bcphi,"dbc2bcphi/F");
  
  
  tree->Branch("bce2nd",&bce2nd,"bce2nd/F");
  tree->Branch("dbceta",&dbceta,"dbceta/F");
  tree->Branch("dbcphi",&dbcphi,"dbcphi/F");
  tree->Branch("bce",&bce,"bce/F");
  tree->Branch("bc2e",&bc2e,"bc2e/F");
  tree->Branch("bc2emax",&bc2emax,"bc2emax/F");
  tree->Branch("bc2e2nd",&bc2e2nd,"bc2e2nd/F");
  tree->Branch("bc2eleft",&bc2eleft,"bc2eleft/F");
  tree->Branch("bc2eright",&bc2eright,"bc2eright/F");
  tree->Branch("bc2etop",&bc2etop,"bc2etop/F");
  tree->Branch("bc2ebottom",&bc2ebottom,"bc2ebottom/F");
  tree->Branch("dbc2eta",&dbc2eta,"dbc2eta/F");
  tree->Branch("dbc2phi",&dbc2phi,"dbc2phi/F");
  tree->Branch("bc2sigietaieta",&bc2sigietaieta,"bc2sigietaieta/F");
  tree->Branch("bc2sigietaiphi",&bc2sigietaiphi,"bc2sigietaiphi/F");
  tree->Branch("bc2sigiphiiphi",&bc2sigiphiiphi,"bc2sigiphiiphi/F");
  tree->Branch("bc2e3x3",&bc2e3x3,"bc2e3x3/F");
  tree->Branch("bc2e5x5",&bc2e5x5,"bc2e5x5/F");
  tree->Branch("bclaste",&bclaste,"bclaste/F");
  tree->Branch("bclaste3x3",&bclaste3x3,"bclaste3x3/F");
  tree->Branch("bclaste5x5",&bclaste5x5,"bclaste5x5/F");
  tree->Branch("bclastsigietaieta",&bclastsigietaieta,"bclastsigietaieta/F");
  tree->Branch("bclastsigiphiiphi",&bclastsigiphiiphi,"bclastsigiphiiphi/F");
  tree->Branch("bclastsigietaiphi",&bclastsigietaiphi,"bclastsigietaiphi/F");
  tree->Branch("dbclasteta",&dbclasteta,"dbclasteta/F");
  tree->Branch("dbclastphi",&dbclastphi,"dbclastphi/F");
  tree->Branch("bclast2e",&bclast2e,"bclast2e/F");
  tree->Branch("bclast2e3x3",&bclast2e3x3,"bclast2e3x3/F");
  tree->Branch("bclast2e5x5",&bclast2e5x5,"bclast2e5x5/F");
  tree->Branch("bclast2sigietaieta",&bclast2sigietaieta,"bclast2sigietaieta/F");
  tree->Branch("bclast2sigiphiiphi",&bclast2sigiphiiphi,"bclast2sigiphiiphi/F");
  tree->Branch("bclast2sigietaiphi",&bclast2sigietaiphi,"bclast2sigietaiphi/F");
  tree->Branch("dbclast2eta",&dbclast2eta,"dbclast2eta/F");
  tree->Branch("dbclast2phi",&dbclast2phi,"dbclast2phi/F");
  tree->Branch("bcseedetacry",&bcseedetacry,"bcseedetacry/F");
  tree->Branch("bc2etacry",&bc2etacry,"bc2etacry/F");
  tree->Branch("bclast2etacry",&bclast2etacry,"bclast2etacry/F");
  tree->Branch("bclastetacry",&bclastetacry,"bclastetacry/F");
  tree->Branch("bcseedphicry",&bcseedphicry,"bcseedphicry/F");
  tree->Branch("bc2phicry",&bc2phicry,"bc2phicry/F");
  tree->Branch("bclast2phicry",&bclast2phicry,"bclast2phicry/F");
  tree->Branch("bclastphicry",&bclastphicry,"bclastphicry/F");
  tree->Branch("bc2ieta",&bc2ieta,"bc2ieta/I");
  tree->Branch("bc2iphi",&bc2iphi,"bc2iphi/I");

}
