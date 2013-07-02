#include "HggEGEnergyCorrector.hh"
#include "HggPhysUtils.cc"
#include "ReadConfig.hh"
#include <string>
#include <iostream>

#define debugEGEnergy 0
HggEGEnergyCorrector::HggEGEnergyCorrector(VecbosBase *r,string cfgFile,Bool_t realData):
  usePhoton(true),
  version(""),
  base(r)
{
  isRealData = realData;
  configFile = cfgFile;
  this->Init();
}

void HggEGEnergyCorrector::getPhotonEnergyCorrection(VecbosPho& pho, bool rescale){
  std::pair<double,double> corr(pho.energy,0);
  if(version.compare("May2012")==0) corr= this->photonEnergyCorrector_May2012(pho,rescale);
  if(version.compare("2011")==0) corr= this->photonEnergyCorrector_CorrectedEnergyWithErrorv2(pho);
  
  pho.correctedEnergy = corr.first;
  pho.correctedEnergyError = corr.second;
}
void HggEGEnergyCorrector::getElectronEnergyCorrection(VecbosEle& ele, bool rescale){
  std::pair<double,double> corr(ele.energy,0);

  if(version.compare("May2012")==0) corr = this->electronEnergyCorrector_May2012(ele,rescale);
  if(version.compare("2011")==0) corr = this->electronEnergyCorrector_CorrectedEnergyWithErrorv2(ele);
  
  ele.correctedEnergy = corr.first;
  ele.correctedEnergyError = corr.second;
}


void HggEGEnergyCorrector::Init(){
  //  fVals = new Float_t[18];
  fVals = new Float_t[73];

  ReadConfig cfg(configFile);
  version = cfg.getParameter("EnergyCorrectionVersion");
  std::string regweightsPho = cfg.getParameter("EnergyCorrectionWeights");
  std::string regweightsEle = cfg.getParameter("EleEnergyCorrectionWeights");

  std::string regweights = (usePhoton ? regweightsPho : regweightsEle);

  std::cout << "Using Correction Weights: " << regweights << std::endl
	    << "Using Correction Version: " << version << std::endl;
  
  TFile *fgbr = new TFile(regweights.c_str(),"READ");
  fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
  fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");
  
  fReaderee = (GBRForest*)fgbr->Get("EECorrection");
  fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");
  fgbr->Close();

  //load the ECAL geometry information
  ecalGeometry = loadecalGapCoordinates(isRealData);
}

std::pair<double,double> HggEGEnergyCorrector::photonEnergyCorrector_May2012(VecbosPho &pho, bool rescale){
  if(base){
    rho = base->rhoJetsFastJet;
    nPV = base->nPV;
  }
  int index=0;
  fVals[0] = pho.SC.rawE;
  fVals[1] = pho.SC.eta;
  fVals[2] = pho.SC.phi;
  fVals[3] = pho.SC.r9;
  fVals[4] = pho.SC.e5x5/pho.SC.rawE;
  fVals[5] = pho.SC.etaWidth;
  fVals[6] = pho.SC.phiWidth;
  fVals[7] = pho.SC.nBCs;
  fVals[8] = pho.HTowOverE; //H/E tower
  fVals[9] = rho;
  fVals[10] = nPV;

  //seed cluster variables
  VecbosBC BC = pho.SC.BCSeed;
  fVals[11] = BC.eta - pho.SC.eta;
  fVals[12] = DeltaPhi(pho.SC.phi,BC.phi);
  fVals[13] = BC.energy/pho.SC.rawE;
  fVals[14] = pho.SC.e3x3/BC.energy;
  fVals[15] = pho.SC.e5x5/BC.energy;

  fVals[16] = pho.SC.sigmaIEtaIEta;
  fVals[17] = pho.SC.sigmaIPhiIPhi;
  fVals[18] = pho.SC.sigmaIEtaIPhi;

  fVals[19] = pho.SC.eMax/BC.energy;
  fVals[20] = pho.SC.e2nd/BC.energy;
  fVals[21] = pho.SC.eTop/BC.energy;
  fVals[22] = pho.SC.eBottom/BC.energy;
  fVals[23] = pho.SC.eLeft/BC.energy;
  fVals[24] = pho.SC.eRight/BC.energy;

  fVals[25] = pho.SC.e2x5Max/BC.energy;
  fVals[26] = pho.SC.e2x5Top/BC.energy;
  fVals[27] = pho.SC.e2x5Bottom/BC.energy;
  fVals[28] = pho.SC.e2x5Left/BC.energy;
  fVals[29] = pho.SC.e2x5Right/BC.energy;
  
  if(pho.isBarrel()){
    fVals[30] = BC.iEta;
    fVals[31] = BC.iPhi;
    fVals[32] = BC.iEta%5;
    fVals[33] = BC.iPhi%2;
    if( abs(BC.iEta <= 25) ) fVals[34] = BC.iEta % 25;
    else fVals[34] = (BC.iEta - 25*abs(BC.iEta)/BC.iEta) % 20;
    fVals[35] = BC.iPhi%20;
    fVals[36] = BC.etaCrystal;
    fVals[37] = BC.phiCrystal;
  }
  else{
    fVals[30] = pho.SC.esEnergy/pho.SC.rawE;
  }
  //finished filling array

  const Double_t varscale = 1.;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (pho.isBarrel()) {
    den = pho.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = pho.SC.rawE + pho.SC.esEnergy;
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;

  if(rescale){
    if(pho.isBarrel()){
      fVals[3] = 1.0045*pho.SC.r9 + 0.001;
      fVals[5] = 1.04301*pho.SC.etaWidth + 0.000618;
      fVals[6] = 1.00002*pho.SC.phiWidth + 0.000371;
      fVals[14] = fVals[3]*pho.SC.e3x3/BC.energy;
      if(fVals[15] <=1.0)
	fVals[15] = TMath::Min(1.0,1.0022*pho.SC.e5x5/BC.energy);
      
      fVals[4] = fVals[15]*fVals[4];

      fVals[16] = 0.891832*pho.SC.sigmaIEtaIEta + 0.0009133;
      fVals[17] = 0.993*pho.SC.sigmaIPhiIPhi;
      
      fVals[19]*=1.012;
      //fVals[20] = 
      fVals[21]*=0.94;
      fVals[22]*=0.94;
      fVals[23]*=0.94;
      fVals[24]*=0.94;
      fVals[25]*=1.006;
      fVals[26]*=1.09;
      fVals[27]*=1.09;
      fVals[28]*=1.09;
      fVals[29]*=1.09;
    }else{ //endcap
      fVals[3] = 1.0086*fVals[3]- 0.0007;
      fVals[5] = 1.04302*fVals[5] - 0.000618;
      fVals[6] = 1.00002*fVals[6] - 0.000371;
      fVals[14]= fVals[3]*fVals[14];
      if(fVals[15] <=1.0)
	fVals[15] = TMath::Min(1.0,1.0022*fVals[15]);
      
      fVals[4] = fVals[15]*fVals[4];
      fVals[16]*= 0.9947;
      
      fVals[19]*=1.005;
      fVals[20]*=1.02;
      fVals[21]*=0.96;
      fVals[22]*=0.96;
      fVals[23]*=0.96;
      fVals[24]*=0.96;
      fVals[25]*=1.0075;
      fVals[26]*=1.13;
      fVals[27]*=1.13;
      fVals[28]*=1.13;
      fVals[29]*=1.13;
    }
  }

  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}

std::pair<double,double> HggEGEnergyCorrector::electronEnergyCorrector_May2012(VecbosEle &ele,bool rescale){
  if(base){
    rho = base->rhoJetsFastJet;
    nPV = base->nPV;
  }
  int index=0;
  fVals[0] = ele.SC.rawE;
  fVals[1] = ele.SC.eta;
  fVals[2] = ele.SC.phi;
  fVals[3] = ele.SC.r9;
  fVals[4] = ele.SC.e5x5/ele.SC.rawE;
  fVals[5] = ele.SC.etaWidth;
  fVals[6] = ele.SC.phiWidth;
  fVals[7] = ele.SC.nBCs;
  fVals[8] = ele.HoverE; //H/E 
  fVals[9] = rho;
  fVals[10] = nPV;

  //seed cluster variables
  VecbosBC BC = ele.SC.BCSeed;
  fVals[11] = BC.eta - ele.SC.eta;
  fVals[12] = DeltaPhi(ele.SC.phi,BC.phi);
  fVals[13] = BC.energy/ele.SC.rawE;
  fVals[14] = ele.SC.e3x3/BC.energy;
  fVals[15] = ele.SC.e5x5/BC.energy;

  fVals[16] = ele.SC.sigmaIEtaIEta;
  fVals[17] = ele.SC.sigmaIPhiIPhi;
  fVals[18] = ele.SC.sigmaIEtaIPhi;

  fVals[19] = ele.SC.eMax/BC.energy;
  fVals[20] = ele.SC.e2nd/BC.energy;
  fVals[21] = ele.SC.eTop/BC.energy;
  fVals[22] = ele.SC.eBottom/BC.energy;
  fVals[23] = ele.SC.eLeft/BC.energy;
  fVals[24] = ele.SC.eRight/BC.energy;

  fVals[25] = ele.SC.e2x5Max/BC.energy;
  fVals[26] = ele.SC.e2x5Top/BC.energy;
  fVals[27] = ele.SC.e2x5Bottom/BC.energy;
  fVals[28] = ele.SC.e2x5Left/BC.energy;
  fVals[29] = ele.SC.e2x5Right/BC.energy;
  
  if( fabs(ele.SC.eta)<1.48 ){
    fVals[30] = BC.iEta;
    fVals[31] = BC.iPhi;
    fVals[32] = BC.iEta%5;
    fVals[33] = BC.iPhi%2;
    if( abs(BC.iEta <= 25) ) fVals[34] = BC.iEta % 25;
    else fVals[34] = (BC.iEta - 25*abs(BC.iEta)/BC.iEta) % 20;
    fVals[35] = BC.iPhi%20;
    fVals[36] = BC.etaCrystal;
    fVals[37] = BC.phiCrystal;
  }
  else{
    fVals[30] = ele.SC.esEnergy/ele.SC.rawE;
  }
  //finished filling array

  const Double_t varscale = 1.;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if( fabs(ele.SC.eta)<1.48 ) {
    den = ele.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = ele.SC.rawE + ele.SC.esEnergy;
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;

  if(rescale){
    if( fabs(ele.SC.eta) < 1.48 ){
      fVals[3] = 1.0045*ele.SC.r9 + 0.001;
      fVals[5] = 1.04301*ele.SC.etaWidth + 0.000618;
      fVals[6] = 1.00002*ele.SC.phiWidth + 0.000371;
      fVals[14] = fVals[3]*ele.SC.e3x3/BC.energy;
      if(fVals[15] <=1.0)
	fVals[15] = TMath::Min(1.0,1.0022*ele.SC.e5x5/BC.energy);
      
      fVals[4] = fVals[15]*fVals[4];

      fVals[16] = 0.891832*ele.SC.sigmaIEtaIEta + 0.0009133;
      fVals[17] = 0.993*ele.SC.sigmaIPhiIPhi;
      
      fVals[19]*=1.012;
      //fVals[20] = 
      fVals[21]*=0.94;
      fVals[22]*=0.94;
      fVals[23]*=0.94;
      fVals[24]*=0.94;
      fVals[25]*=1.006;
      fVals[26]*=1.09;
      fVals[27]*=1.09;
      fVals[28]*=1.09;
      fVals[29]*=1.09;
    }else{ //endcap
      fVals[3] = 1.0086*fVals[3]- 0.0007;
      fVals[5] = 1.04302*fVals[5] - 0.000618;
      fVals[6] = 1.00002*fVals[6] - 0.000371;
      fVals[14]= fVals[3]*fVals[14];
      if(fVals[15] <=1.0)
	fVals[15] = TMath::Min(1.0,1.0022*fVals[15]);
      
      fVals[4] = fVals[15]*fVals[4];
      fVals[16]*= 0.9947;
      
      fVals[19]*=1.005;
      fVals[20]*=1.02;
      fVals[21]*=0.96;
      fVals[22]*=0.96;
      fVals[23]*=0.96;
      fVals[24]*=0.96;
      fVals[25]*=1.0075;
      fVals[26]*=1.13;
      fVals[27]*=1.13;
      fVals[28]*=1.13;
      fVals[29]*=1.13;
    }
  }

  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}

std::pair<double,double> HggEGEnergyCorrector::CorrectedEnergyWithError(int j){
  VecbosPho pho(base,j); // initialize the photon with index j

  int iSC = base->superClusterIndexPho[j]; // supercluster index for this photon

  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry,pho.SC.eta,pho.SC.phi);
  
  fVals[0]  = pho.SC.rawE;
  fVals[1]  = pho.SC.r9;
  fVals[2]  = pho.SC.eta;
  fVals[3]  = pho.SC.phi;
  fVals[4]  = pho.SC.e5x5/pho.SC.rawE;
  
  Bool_t isbarrel = fabs(pho.SC.eta) < 1.48; 
  
  if( isbarrel){
    fVals[5]  = thisGeometry._aC; 
    fVals[6]  = thisGeometry._aS;
    fVals[7]  = thisGeometry._aM;
    fVals[8]  = thisGeometry._bC;
    fVals[9]  = thisGeometry._bS;
    fVals[10] = thisGeometry._bM;
    fVals[11] = pho.HoverE;
    fVals[12] = pho.SC.etaWidth;
    fVals[13] = pho.SC.phiWidth;
    fVals[14] = pho.SC.sigmaIEtaIEta;
  }else{
    fVals[5]  = pho.SC.rawE; //photonscpreshowerEnergy[j] / photonscrawEnergy[j]; /// DO WE STORE THIS????
    fVals[6]  = thisGeometry._xZ;
    fVals[7]  = thisGeometry._aC;
    fVals[8]  = thisGeometry._aS;
    fVals[9]  = thisGeometry._aM;
    fVals[10] = thisGeometry._yZ;
    fVals[11] = thisGeometry._bC;
    fVals[12] = thisGeometry._bS;
    fVals[13] = thisGeometry._bM;
    fVals[14] = pho.SC.HoverE;
    fVals[15] = pho.SC.etaWidth;
    fVals[16] = pho.SC.phiWidth;
    fVals[17] = pho.SC.sigmaIEtaIEta;
  }

  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = pho.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = pho.SC.rawE + pho.SC.esEnergy;
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}



std::pair<double,double> HggEGEnergyCorrector::electronEnergyCorrector_CorrectedEnergyWithError(int j){
  VecbosEle ele(base,j);
    
  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry,ele.SC.eta,ele.SC.phi);
  
  fVals[0]  = ele.SC.rawE;
  fVals[1]  = ele.SC.r9;
  fVals[2]  = ele.SC.eta;
  fVals[3]  = ele.SC.phi;
  fVals[4]  = ele.SC.e5x5/ele.SC.rawE;
  
  Bool_t isbarrel = fabs(ele.SC.eta) < 1.48; 
  
  if( isbarrel){
    fVals[5]  = thisGeometry._aC; 
    fVals[6]  = thisGeometry._aS;
    fVals[7]  = thisGeometry._aM;
    fVals[8]  = thisGeometry._bC;
    fVals[9]  = thisGeometry._bS;
    fVals[10] = thisGeometry._bM;
    fVals[11] = ele.SC.HoverE;
    fVals[12] = ele.SC.etaWidth;
    fVals[13] = ele.SC.phiWidth;
    fVals[14] = ele.SC.sigmaIEtaIEta;

  }else{
    fVals[5]  = ele.esEnergy/ele.SC.rawE;
    fVals[6]  = thisGeometry._xZ;
    fVals[7]  = thisGeometry._aC;
    fVals[8]  = thisGeometry._aS;
    fVals[9]  = thisGeometry._aM;
    fVals[10] = thisGeometry._yZ;
    fVals[11] = thisGeometry._bC;
    fVals[12] = thisGeometry._bS;
    fVals[13] = thisGeometry._bM;
    fVals[14] = ele.SC.HoverE;
    fVals[15] = ele.SC.etaWidth;
    fVals[16] = ele.SC.phiWidth;
    fVals[17] = ele.SC.sigmaIEtaIEta;
  }
  
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = ele.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = ele.SC.rawE + ele.esEnergy;
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}


std::pair<double,double> HggEGEnergyCorrector::electronEnergyCorrector_CorrectedEnergyWithErrorv2(VecbosEle &ele){
    return std::pair<double,double>(ele.SC.energy,0);
#if false 
  ///no correction for electron of tracker-deriven seed only ( not used in analysis anyway)
  if( !ele.isEcalDriven ){
    return std::pair<double,double>(ele.SC.energy,0);
  }

  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry, ele.SC.eta, ele.SC.phi);
  
 

  Bool_t isbarrel = fabs(ele.SC.eta) < 1.48; 
  
  if( isbarrel){
    

    fVals[0]  = ele.SC.rawE;
    fVals[1]  = ele.SC.r9;
    fVals[2]  = ele.SC.eta;
    fVals[3]  = ele.SC.phi;
    fVals[4]  = ele.SC.e5x5/ele.SC.rawE;
    fVals[5] = ele.HoverE;
    fVals[6] = ele.SC.etaWidth;
    fVals[7] = ele.SC.phiWidth;
    ///bc
    VecbosBC BC1 = ele.SC.basicClusters[0];
    
    fVals[8]  = BC1.eta - ele.SC.eta;
    fVals[9]  = DeltaPhi(BC1.phi,ele.SC.phi);
    fVals[10]  = BC1.energy/ele.SC.rawE;
    fVals[11]  = BC1.e3x3/BC1.energy;
    fVals[12]  = BC1.e5x5/BC1.energy;
    fVals[13] = BC1.sigmaIEtaIEta;
    fVals[14] = BC1.sigmaIPhiIPhi;
    fVals[15] = BC1.sigmaIEtaIPhi;
    fVals[16] = BC1.eMax/BC1.energy;
    fVals[17] = log(BC1.e2nd/BC1.eMax);
    fVals[18] = log(BC1.eTop/BC1.eMax);
    fVals[19] = log(BC1.eBottom/BC1.eMax);
    fVals[20] = log(BC1.eLeft/BC1.eMax);
    fVals[21] = log(BC1.eRight/BC1.eMax);
    fVals[22] = (BC1.eTop-BC1.eBottom)/(BC1.eTop + BC1.eBottom);
    fVals[23] = (BC1.eLeft-BC1.eRight)/(BC1.eLeft+BC1.eRight);

    if(ele.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = ele.SC.basicClusters[1];
      fVals[24]  = BC2.eta - ele.SC.eta;
      fVals[25]  = DeltaPhi(BC2.phi,ele.SC.phi);
      fVals[26]  = BC2.energy/ele.SC.rawE;
      fVals[27]  = BC2.e3x3/BC2.energy;
      fVals[28]  = BC2.e5x5/BC2.energy;
      fVals[29] = BC2.sigmaIEtaIEta;
      fVals[30] = BC2.sigmaIPhiIPhi;
      fVals[31] = BC2.sigmaIEtaIPhi;
      fVals[32] = BC2.eMax/BC2.energy;
      fVals[33] = log(BC2.e2nd/BC2.eMax);
      fVals[34] = log(BC2.eTop/BC2.eMax);
      fVals[35] = log(BC2.eBottom/BC2.eMax);
      fVals[36] = log(BC2.eLeft/BC2.eMax);
      fVals[37] = log(BC2.eRight/BC2.eMax);
      fVals[38] = (BC2.eTop-BC2.eBottom)/(BC2.eTop + BC2.eBottom);
      fVals[39] = (BC2.eLeft-BC2.eRight)/(BC2.eLeft+BC2.eRight);
    }else{ // only one basic cluster
      for(int i=24;i<40;i++) fVals[i] = 0.;
    }
    
    if(ele.SC.basicClusters.size() > 2){ // now look for the lowest energy basic cluster (for pileup mitigation)
      VecbosBC BCL = ele.SC.basicClusters.back(); // last element in the energy sorted list
      fVals[40] = BCL.eta-ele.SC.eta;
      fVals[41] = DeltaPhi(BCL.phi,ele.SC.phi);
      fVals[42] = BCL.energy/ele.SC.rawE;
      fVals[43] = BCL.e3x3/BCL.energy;
      fVals[44] = BCL.e5x5/BCL.energy;
      fVals[45] = BCL.sigmaIEtaIEta;
      fVals[46] = BCL.sigmaIPhiIPhi;
      fVals[47] = BCL.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=40;i<48;i++) fVals[i] = 0;
    }
    
    if(ele.SC.basicClusters.size() > 3){ // now look for the 2nd lowest energy basic cluster (for pileup mitigation)
      VecbosBC BC2L = ele.SC.basicClusters.at(ele.SC.basicClusters.size()-2); // last element in the energy sorted list
      fVals[48] = BC2L.eta-ele.SC.eta;
      fVals[49] = DeltaPhi(BC2L.phi,ele.SC.phi);
      fVals[50] = BC2L.energy/ele.SC.rawE;
      fVals[51] = BC2L.e3x3/BC2L.energy;
      fVals[52] = BC2L.e5x5/BC2L.energy;
      fVals[53] = BC2L.sigmaIEtaIEta;
      fVals[54] = BC2L.sigmaIPhiIPhi;
      fVals[55] = BC2L.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=48;i<56;i++) fVals[i] = 0;
    }
    
   
 //local coordinates and crystal indices
    //seed cluster

    fVals[56] = BC1.iEta; //crystal ieta
    fVals[57] = BC1.iPhi; //crystal iphi
    fVals[58] = BC1.iEta%5; //submodule boundary eta symmetry
    fVals[59] = BC1.iPhi%2; //submodule boundary phi symmetry    
    fVals[60] = (TMath::Abs(BC1.iEta)<=25)*(BC1.iEta%25) + (TMath::Abs(BC1.iEta)>25)*((BC1.iEta-25*TMath::Abs(BC1.iEta)/BC1.iEta)%20);  //module boundary eta approximate symmetry            
    fVals[61] = BC1.iPhi%20; //module boundary phi symmetry                                                                                                                    
    fVals[62] = BC1.etaCrystal; //local coordinates with respect to closest crystal center at nominal shower depth                                                                 
    fVals[63] = BC1.phiCrystal;

    //2nd cluster (meaningful gap corrections for converted photons)                                                                                                        
    if(ele.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = ele.SC.basicClusters[1];
      fVals[64] = BC2.iEta;
      fVals[65] = BC2.iPhi;
      fVals[66] = BC2.iEta%5;
      fVals[67] = BC2.iPhi%2;
      fVals[68] = (TMath::Abs(BC2.iEta)<=25)*(BC2.iEta%25) + (TMath::Abs(BC2.iEta)>25)*((BC2.iEta-25*TMath::Abs(BC2.iEta)/BC2.iEta)%20);
      fVals[69] = BC2.iPhi%20;
      fVals[70] = BC2.etaCrystal;
      fVals[71] = BC2.phiCrystal;
    }else{
      for(int i=64;i<72;i++) fVals[i] =0.;
    }
    //Nb. of vertex
    fVals[72] = base->nPV;
    
    
    
  }else{
    
    fVals[0]  = ele.SC.rawE;
    fVals[1]  = ele.SC.r9;
    fVals[2]  = ele.SC.eta;
    fVals[3]  = ele.SC.phi;
    fVals[4]  = ele.SC.e5x5/ele.SC.rawE;
    fVals[5]  = ele.SC.etaWidth;
    fVals[6]  = ele.SC.phiWidth;
    //Nb. of vertex
    fVals[7] = base->nPV; 
    
  }
  
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = ele.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = ele.SC.rawE + ele.esEnergy;
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  
  return std::pair<double,double>(ecor,ecorerr);
  
#endif  
}




std::pair<double,double> HggEGEnergyCorrector::photonEnergyCorrector_CorrectedEnergyWithErrorv2(VecbosPho &pho){
  return std::pair<double,double>(0,0);
#if false  // new vecbos version totally breaks compatibility
  if(debugEGEnergy) cout << "Getting ECAL Coordinates" << endl;
  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry, pho.SC.eta, pho.SC.phi);
  if(debugEGEnergy) cout << "Done" << endl;
  
 

  Bool_t isbarrel = fabs(pho.SC.eta) < 1.48; 
  
  if( isbarrel){
        
    fVals[0]  = pho.SC.rawE;
    fVals[1]  = pho.SC.r9;
    fVals[2]  = pho.SC.eta;
    fVals[3]  = pho.SC.phi;
    fVals[4]  = pho.SC.e5x5/pho.SC.rawE;
    fVals[5] = pho.HoverE;
    fVals[6] = pho.SC.etaWidth;
    fVals[7] = pho.SC.phiWidth;
    ///bc
    VecbosBC BC1 = pho.SC.basicClusters[0];
    if(debugEGEnergy) cout << "Basic Cluster 1" <<endl;
    if(debugEGEnergy){
      cout << BC1.eta - pho.SC.eta << endl
	   << DeltaPhi(BC1.phi,pho.SC.phi) << endl
	   << BC1.energy/pho.SC.rawE << endl
	   << BC1.e3x3/BC1.energy << endl
	   << BC1.e5x5/BC1.energy << endl
	   << sqrt(BC1.sigmaIEtaIEta) << endl
	   << sqrt(BC1.sigmaIPhiIPhi) << endl
	   << BC1.sigmaIEtaIPhi << endl 
	   << log(BC1.e2nd/BC1.eMax) << endl
	   <<log(BC1.eTop/BC1.eMax) << endl
	   <<log(BC1.eLeft/BC1.eMax) << endl
	   << log(BC1.eRight/BC1.eMax) << endl
	   << (BC1.eTop-BC1.eBottom)/(BC1.eTop + BC1.eBottom) << endl
	   << (BC1.eLeft-BC1.eRight)/(BC1.eLeft+BC1.eRight) << endl;
    }
    fVals[8]  = BC1.eta - pho.SC.eta;
    fVals[9]  = DeltaPhi(BC1.phi,pho.SC.phi);
    fVals[10]  = BC1.energy/pho.SC.rawE;
    fVals[11]  = BC1.e3x3/BC1.energy;
    fVals[12]  = BC1.e5x5/BC1.energy;
    fVals[13] = BC1.sigmaIEtaIEta;
    fVals[14] = BC1.sigmaIPhiIPhi;
    fVals[15] = BC1.sigmaIEtaIPhi;
    fVals[16] = BC1.eMax/BC1.energy;
    fVals[17] = log(BC1.e2nd/BC1.eMax);
    fVals[18] = log(BC1.eTop/BC1.eMax);
    fVals[19] = log(BC1.eBottom/BC1.eMax);
    fVals[20] = log(BC1.eLeft/BC1.eMax);
    fVals[21] = log(BC1.eRight/BC1.eMax);
    fVals[22] = (BC1.eTop-BC1.eBottom)/(BC1.eTop + BC1.eBottom);
    fVals[23] = (BC1.eLeft-BC1.eRight)/(BC1.eLeft+BC1.eRight);

    if(pho.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = pho.SC.basicClusters[1];
      if(debugEGEnergy) cout << "Basic Cluster 2" <<endl;
    if(debugEGEnergy){
      cout << BC2.eta - pho.SC.eta << endl
	   << DeltaPhi(BC2.phi,pho.SC.phi) << endl
	   << BC2.energy/pho.SC.rawE << endl
	   << BC2.e3x3/BC2.energy << endl
	   << BC2.e5x5/BC2.energy << endl
	   << sqrt(BC2.sigmaIEtaIEta) << endl
	   << sqrt(BC2.sigmaIPhiIPhi) << endl
	   << BC2.sigmaIEtaIPhi << endl 
	   << log(BC2.e2nd/BC2.eMax) << endl
	   <<log(BC2.eTop/BC2.eMax) << endl
	   <<log(BC2.eLeft/BC2.eMax) << endl
	   << log(BC2.eRight/BC2.eMax) << endl
	   << (BC2.eTop-BC2.eBottom)/(BC2.eTop + BC2.eBottom) << endl
	   << (BC2.eLeft-BC2.eRight)/(BC2.eLeft+BC2.eRight) << endl;
    }
      fVals[24]  = BC2.eta - pho.SC.eta;
      fVals[25]  = DeltaPhi(BC2.phi,pho.SC.phi);
      fVals[26]  = BC2.energy/pho.SC.rawE;
      fVals[27]  = BC2.e3x3/BC2.energy;
      fVals[28]  = BC2.e5x5/BC2.energy;
      fVals[29] = BC2.sigmaIEtaIEta;
      fVals[30] = BC2.sigmaIPhiIPhi;
      fVals[31] = BC2.sigmaIEtaIPhi;
      fVals[32] = BC2.eMax/BC2.energy;
      fVals[33] = log(BC2.e2nd/BC2.eMax);
      fVals[34] = log(BC2.eTop/BC2.eMax);
      fVals[35] = log(BC2.eBottom/BC2.eMax);
      fVals[36] = log(BC2.eLeft/BC2.eMax);
      fVals[37] = log(BC2.eRight/BC2.eMax);
      fVals[38] = (BC2.eTop-BC2.eBottom)/(BC2.eTop + BC2.eBottom);
      fVals[39] = (BC2.eLeft-BC2.eRight)/(BC2.eLeft+BC2.eRight);
    }else{ // only one basic cluster
      for(int i=24;i<40;i++) fVals[i] = 0.;
    }
    
    if(pho.SC.basicClusters.size() > 2){ // now look for the lowest energy basic cluster (for pileup mitigation)
      VecbosBC BCL = pho.SC.basicClusters.back(); // last element in the energy sorted list
      if(debugEGEnergy) cout << "Basic Cluster L" <<endl;
    if(debugEGEnergy){
      cout << BCL.eta - pho.SC.eta << endl
	   << DeltaPhi(BCL.phi,pho.SC.phi) << endl
	   << BCL.energy/pho.SC.rawE << endl
	   << BCL.e3x3/BCL.energy << endl
	   << BCL.e5x5/BCL.energy << endl
	   << sqrt(BCL.sigmaIEtaIEta) << endl
	   << sqrt(BCL.sigmaIPhiIPhi) << endl
	   << BCL.sigmaIEtaIPhi << endl 
	   << log(BCL.e2nd/BCL.eMax) << endl
	   <<log(BCL.eTop/BCL.eMax) << endl
	   <<log(BCL.eLeft/BCL.eMax) << endl
	   << log(BCL.eRight/BCL.eMax) << endl
	   << (BCL.eTop-BCL.eBottom)/(BCL.eTop + BCL.eBottom) << endl
	   << (BCL.eLeft-BCL.eRight)/(BCL.eLeft+BCL.eRight) << endl;
    }
      fVals[40] = BCL.eta-pho.SC.eta;
      fVals[41] = DeltaPhi(BCL.phi,pho.SC.phi);
      fVals[42] = BCL.energy/pho.SC.rawE;
      fVals[43] = BCL.e3x3/BCL.energy;
      fVals[44] = BCL.e5x5/BCL.energy;
      fVals[45] = BCL.sigmaIEtaIEta;
      fVals[46] = BCL.sigmaIPhiIPhi;
      fVals[47] = BCL.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=40;i<48;i++) fVals[i] = 0;
    }
    
    if(pho.SC.basicClusters.size() > 3){ // now look for the 2nd lowest energy basic cluster (for pileup mitigation)
      VecbosBC BC2L = pho.SC.basicClusters.at(pho.SC.basicClusters.size()-2); // last element in the energy sorted list
      if(debugEGEnergy) cout << "Basic Cluster PU" <<endl;
    if(debugEGEnergy){
      cout << BC2L.eta - pho.SC.eta << endl
	   << DeltaPhi(BC2L.phi,pho.SC.phi) << endl
	   << BC2L.energy/pho.SC.rawE << endl
	   << BC2L.e3x3/BC2L.energy << endl
	   << BC2L.e5x5/BC2L.energy << endl
	   << sqrt(BC2L.sigmaIEtaIEta) << endl
	   << sqrt(BC2L.sigmaIPhiIPhi) << endl
	   << BC2L.sigmaIEtaIPhi << endl 
	   << log(BC2L.e2nd/BC2L.eMax) << endl
	   <<log(BC2L.eTop/BC2L.eMax) << endl
	   <<log(BC2L.eLeft/BC2L.eMax) << endl
	   << log(BC2L.eRight/BC2L.eMax) << endl
	   << (BC2L.eTop-BC2L.eBottom)/(BC2L.eTop + BC2L.eBottom) << endl
	   << (BC2L.eLeft-BC2L.eRight)/(BC2L.eLeft+BC2L.eRight) << endl;
    }
      fVals[48] = BC2L.eta-pho.SC.eta;
      fVals[49] = DeltaPhi(BC2L.phi,pho.SC.phi);
      fVals[50] = BC2L.energy/pho.SC.rawE;
      fVals[51] = BC2L.e3x3/BC2L.energy;
      fVals[52] = BC2L.e5x5/BC2L.energy;
      fVals[53] = BC2L.sigmaIEtaIEta;
      fVals[54] = BC2L.sigmaIPhiIPhi;
      fVals[55] = BC2L.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=48;i<56;i++) fVals[i] = 0;
    }
    if(debugEGEnergy) cout << "Basic Cluster 1 iEtas" <<endl;
    fVals[56] = BC1.iEta; //crystal ieta
    fVals[57] = BC1.iPhi; //crystal iphi
    if(debugEGEnergy) cout << "modules" <<endl;
    fVals[58] = (BC1.iEta)%5; //submodule boundary eta symmetry
    fVals[59] = (BC1.iPhi)%2; //submodule boundary phi symmetry    
    if(debugEGEnergy) cout << "TMath stuff" <<endl;
    fVals[60] = (TMath::Abs(BC1.iEta)<=25 ? (BC1.iEta)%25 : (BC1.iEta-25*(BC1.iEta>0))%20 );
      // (TMath::Abs(BC1.iEta)<=25)*(BC1.iEta%25) + (TMath::Abs(BC1.iEta)>25)*((BC1.iEta-25*TMath::Abs(BC1.iEta)/BC1.iEta)%20);  //module boundary eta approximate symmetry            
    if(debugEGEnergy) cout << "not much else" <<endl;
    fVals[61] = (BC1.iPhi)%20; //module boundary phi symmetry                                                                                                                    
    fVals[62] = BC1.etaCrystal; //local coordinates with respect to closest crystal center at nominal shower depth                                                                 
    fVals[63] = BC1.phiCrystal;

    //2nd cluster (meaningful gap corrections for converted photons)                                                                                                        
    if(pho.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = pho.SC.basicClusters[1];
      if(debugEGEnergy) cout << "Basic Cluster 2 iEtas" <<endl;
      fVals[64] = BC2.iEta;
      fVals[65] = BC2.iPhi;
      fVals[66] = BC2.iEta%5;
      fVals[67] = BC2.iPhi%2;
      fVals[68] =  (TMath::Abs(BC2.iEta)<=25 ? (BC2.iEta)%25 : (BC2.iEta-25*(BC2.iEta>0))%20 );
	//(TMath::Abs(BC2.iEta)<=25)*(BC2.iEta%25) + (TMath::Abs(BC2.iEta)>25)*((BC2.iEta-25*TMath::Abs(BC2.iEta)/BC2.iEta)%20);
      fVals[69] = BC2.iPhi%20;
      fVals[70] = BC2.etaCrystal;
      fVals[71] = BC2.phiCrystal;
    }else{
      for(int i=64;i<72;i++) fVals[i] =0.;
    }
    //Nb. of vertex
    fVals[72] = base->nPV;
    
    
    
  }else{
    
    fVals[0]  = pho.SC.rawE;
    fVals[1]  = pho.SC.r9;
    fVals[2]  = pho.SC.eta;
    fVals[3]  = pho.SC.phi;
    fVals[4]  = pho.SC.e5x5/pho.SC.rawE;
    fVals[5]  = pho.SC.etaWidth;
    fVals[6]  = pho.SC.phiWidth;
    //Nb. of vertex
    fVals[7] = base->nPV; 
    
  }
  if(debugEGEnergy) cout << "doing GBR stuff" <<endl;

  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = pho.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = pho.SC.rawE + pho.SC.esEnergy; 
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  if(debugEGEnergy) cout << "done returning. ecor: " << ecor/den << "  +-  " << ecorerr/den <<endl;

  return std::pair<double,double>(ecor,ecorerr);
  
#endif  
}
