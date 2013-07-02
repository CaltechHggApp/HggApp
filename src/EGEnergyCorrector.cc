GBRForest *fReadereb;
GBRForest *fReaderebvariance;


GBRForest *fReaderee;
GBRForest *fReadereevariance;

Float_t *fVals;


void EGEnergyCorrector_Init(){
  //  fVals = new Float_t[18];
    fVals = new Float_t[73];
    
  string regweights; 
  if( testPhotonCorr == 99){
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrph.root";
  }else if( testPhotonCorr == 98) {
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrele.root";
  }else if( testPhotonCorr == 97){
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrv2ele.root";
  }else if( testPhotonCorr == 96){
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrv2ph.root";
  }
  

  
  
  cout<<"regweights " << regweights.c_str()<<endl; 
  
  TFile *fgbr = new TFile(regweights.c_str(),"READ");
  fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
  fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");
  
  fReaderee = (GBRForest*)fgbr->Get("EECorrection");
  fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");
  fgbr->Close();
  
    
  cout<<"regweights loaded " <<endl;
  

}



std::pair<double,double> photonEnergyCorrector_CorrectedEnergyWithError(int j){
    
  getGapCoordinates(photonsceta[j],photonscphi[j]);
  
  fVals[0]  = photonscrawEnergy[j];
  fVals[1]  = photonr9[j];
  fVals[2]  = photonsceta[j];
  fVals[3]  = photonscphi[j];
  fVals[4]  = photone5x5[j]/photonscrawEnergy[j];
  
  Bool_t isbarrel = fabs(photonsceta[j]) < 1.48; 
  
  if( isbarrel){
    fVals[5]  = _aC; 
    fVals[6]  = _aS;
    fVals[7]  = _aM;
    fVals[8]  = _bC;
    fVals[9]  = _bS;
    fVals[10] = _bM;
    fVals[11] = photonhadronicOverEm[j];
    fVals[12] = photonscetaWidth[j];
    fVals[13] = photonscphiWidth[j];
    fVals[14] = photonsigmaIetaIeta[j];
  }else{
    fVals[5]  = photonscpreshowerEnergy[j] / photonscrawEnergy[j];
    fVals[6]  = _xZ;
    fVals[7]  = _aC;
    fVals[8]  = _aS;
    fVals[9]  = _aM;
    fVals[10] = _yZ;
    fVals[11] = _bC;
    fVals[12] = _bS;
    fVals[13] = _bM;
    fVals[14] = photonhadronicOverEm[j];
    fVals[15] = photonscetaWidth[j];
    fVals[16] = photonscphiWidth[j];
    fVals[17] = photonsigmaIetaIeta[j];
  }

  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = photonscrawEnergy[j];
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = photonscrawEnergy[j] + photonscpreshowerEnergy[j]; 
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
  
  
  
}



std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithError(int j){
    
  getGapCoordinates(electronsceta[j],electronscphi[j]);
  
  fVals[0]  = electronscrawEnergy[j];
  fVals[1]  = electrone3x3[j]/electronscrawEnergy[j];
  fVals[2]  = electronsceta[j];
  fVals[3]  = electronscphi[j];
  fVals[4]  = electrone5x5[j]/electronscrawEnergy[j];
  
  Bool_t isbarrel = fabs(electronsceta[j]) < 1.48; 
  
  if( isbarrel){
    fVals[5]  = _aC; 
    fVals[6]  = _aS;
    fVals[7]  = _aM;
    fVals[8]  = _bC;
    fVals[9]  = _bS;
    fVals[10] = _bM;
    fVals[11] = electronhcalOverEcal[j];
    fVals[12] = electronscetaWidth[j];
    fVals[13] = electronscphiWidth[j];
    fVals[14] = electronsigmaIetaIeta[j];
    

  }else{
    fVals[5]  = electronscpreshowerEnergy[j] / electronscrawEnergy[j];
    fVals[6]  = _xZ;
    fVals[7]  = _aC;
    fVals[8]  = _aS;
    fVals[9]  = _aM;
    fVals[10] = _yZ;
    fVals[11] = _bC;
    fVals[12] = _bS;
    fVals[13] = _bM;
    fVals[14] = electronhcalOverEcal[j];
    fVals[15] = electronscetaWidth[j];
    fVals[16] = electronscphiWidth[j];
    fVals[17] = electronsigmaIetaIeta[j];
  }
  
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = electronscrawEnergy[j];
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = electronscrawEnergy[j] + electronscpreshowerEnergy[j]; 
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  

//   if( isbarrel){
//     cout<<"checkele " << greader->GetResponse(fVals) <<endl; 
//     for(int j=0;j<15; j++){
//       cout<<" fVals" << fVals[j]<<endl; 
//     }
//   }


  return std::pair<double,double>(ecor,ecorerr);
  
  
  
}


std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithErrorv2(int j){
    
  
  ///no correction for electron of tracker-deriven seed only ( not used in analysis anyway)
  if( !electronecalDrivenSeed[j]){
    return std::pair<double,double>(electronscenergy[j],0);
  }
  

  getGapCoordinates(electronsceta[j],electronscphi[j]);
  
 

  Bool_t isbarrel = fabs(electronsceta[j]) < 1.48; 
  
  if( isbarrel){
    
    
    fVals[0]  = electronscrawEnergy[j];
    fVals[1]  = electrone3x3[j]/electronscrawEnergy[j];
    fVals[2]  = electronsceta[j];
    fVals[3]  = electronscphi[j];
    fVals[4]  = electrone5x5[j]/electronscrawEnergy[j];
    fVals[5] = electronhcalOverEcal[j];
    fVals[6] = electronscetaWidth[j];
    fVals[7] = electronscphiWidth[j];
        
    ///bc
    float bemax = electronbcseedeMax[j];
    float be2nd = electronbcseede2nd[j];
    float betop = electronbcseedeTop[j];
    float bebottom = electronbcseedeBottom[j];
    float beleft = electronbcseedeLeft[j];
    float beright = electronbcseedeRight[j];
    

    fVals[8]  = electronbcseedeta[j] - electronsceta[j]; 
    fVals[9]  = DeltaPhi(electronbcseedphi[j],electronscphi[j]);
    fVals[10]  = electronbcseedenergy[j]/ electronscrawEnergy[j];
    fVals[11]  = electronbcseede3x3[j] / electronbcseedenergy[j];
    fVals[12]  = electronbcseede5x5[j] / electronbcseedenergy[j];
    fVals[13] = sqrt(electronbcseedCovIEtaIEta[j]);
    fVals[14] = sqrt(electronbcseedCovIPhiIPhi[j]);
    fVals[15] = electronbcseedCovIEtaIPhi[j];
    fVals[16] = bemax/electronbcseedenergy[j];
    fVals[17] = log(be2nd/bemax);
    fVals[18] = log(betop/bemax);
    fVals[19] = log(bebottom/bemax);
    fVals[20] = log(beleft/bemax);
    fVals[21] = log(beright/bemax);
    fVals[22] = (betop-bebottom)/(betop+bebottom);
    fVals[23] = (beleft-beright)/(beleft+beright);

    

    ///bc2
    bool hasbc2 = electronbc2eMax[j] > 0; 
    float bc2emax = electronbc2eMax[j];
    float bc2e2nd = electronbc2e2nd[j];
    float bc2etop = electronbc2eTop[j];
    float bc2ebottom = electronbc2eBottom[j];
    float bc2eleft = electronbc2eLeft[j];
    float bc2eright = electronbc2eRight[j];
    
    fVals[24] = hasbc2 ? electronbc2eta[j] - electronsceta[j]: 0; 
    fVals[25] = hasbc2 ? DeltaPhi(electronbc2phi[j],electronscphi[j]): 0; 
    fVals[26] = hasbc2 ? electronbc2energy[j]/ electronscrawEnergy[j]: 0; 
    fVals[27] = hasbc2 ? electronbc2e3x3[j]/electronbc2energy[j]: 0; 
    fVals[28] = hasbc2 ? electronbc2e5x5[j]/electronbc2energy[j]: 0;
    fVals[29] = hasbc2 ? sqrt(electronbc2CovIEtaIEta[j]): 0;
    fVals[30] = hasbc2 ? sqrt(electronbc2CovIPhiIPhi[j]): 0; 
    
    
    fVals[31] = hasbc2 ? electronbc2CovIEtaIPhi[j]: 0; 
    ///bug ? 
    //fVals[31] = hasbc2 ? electronbcseedCovIEtaIPhi[j]: 0; 
    
    fVals[32] = hasbc2 ? bc2emax/electronbc2energy[j] : 0.;
    fVals[33] = hasbc2 ? log(bc2e2nd/bc2emax) : 0.;
    fVals[34] = hasbc2 ? log(bc2etop/bc2emax) : 0.;
    fVals[35] = hasbc2 ? log(bc2ebottom/bc2emax) : 0.;
    fVals[36] = hasbc2 ? log(bc2eleft/bc2emax) : 0.;
    fVals[37] = hasbc2 ? log(bc2eright/bc2emax) : 0.;
    fVals[38] = hasbc2 ? (bc2etop-bc2ebottom)/(bc2etop+bc2ebottom) : 0.;
    fVals[39] = hasbc2 ? (bc2eleft-bc2eright)/(bc2eleft+bc2eright) : 0.;
    
    //bclast
    bool hasbclast = electronbclastenergy[j]> 0; 
    fVals[40] = hasbclast ? (electronbclasteta[j]-electronsceta[j]) : 0.;
    fVals[41] = hasbclast ? DeltaPhi(electronbclastphi[j],electronscphi[j]) : 0.;
    fVals[42] = hasbclast ? electronbclastenergy[j] / electronscrawEnergy[j] : 0.;  
    fVals[43] = hasbclast ? electronbclaste3x3[j] / electronbclastenergy[j] :0.; 
    fVals[44] = hasbclast ? electronbclaste5x5[j] / electronbclastenergy[j] :0.; 
    fVals[45] = hasbclast ? sqrt(electronbclastCovIEtaIEta[j]): 0;
    fVals[46] = hasbclast ? sqrt(electronbclastCovIPhiIPhi[j]): 0; 
    fVals[47] = hasbclast ? electronbclastCovIEtaIPhi[j]: 0; 
    
    ///bclast2
    bool hasbclast2 = electronbclast2energy[j]> 0; 
    fVals[48] = hasbclast2 ? (electronbclast2eta[j]-electronsceta[j]) : 0.;
    fVals[49] = hasbclast2 ? DeltaPhi(electronbclast2phi[j],electronscphi[j]) : 0.;
    fVals[50] = hasbclast2 ? electronbclast2energy[j] / electronscrawEnergy[j] : 0.;  
    fVals[51] = hasbclast2 ? electronbclast2e3x3[j] / electronbclast2energy[j] :0.; 
    fVals[52] = hasbclast2 ? electronbclast2e5x5[j] / electronbclast2energy[j] :0.; 
    fVals[53] = hasbclast2 ? sqrt(electronbclast2CovIEtaIEta[j]): 0;
    fVals[54] = hasbclast2 ? sqrt(electronbclast2CovIPhiIPhi[j]): 0; 
    fVals[55] = hasbclast2 ? electronbclast2CovIEtaIPhi[j]: 0; 
    
    //local coordinates and crystal indices                                                                                                                                 
    //seed cluster                                                                                                                                                          
    float betacry, bphicry, bthetatilt, bphitilt;
    int bieta, biphi;
    
    bieta = electronbcseedieta[j];
    biphi = electronbcseediphi[j];
    betacry = electronbcseedetacry[j];
    bphicry = electronbcseedphicry[j];
    
    fVals[56] = bieta; //crystal ieta                                                                                                                                       
    fVals[57] = biphi; //crystal iphi                                                                                                                                       
    fVals[58] = bieta%5; //submodule boundary eta symmetry                                                                                                                  
    fVals[59] = biphi%2; //submodule boundary phi symmetry    
    fVals[60] = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry            
    fVals[61] = biphi%20; //module boundary phi symmetry                                                                                                                    
    fVals[62] = betacry; //local coordinates with respect to closest crystal center at nominal shower depth                                                                 
    fVals[63] = bphicry;

    //2nd cluster (meaningful gap corrections for converted photons)                                                                                                        
    float bc2etacry, bc2phicry, bc2thetatilt, bc2phitilt;
    int bc2ieta, bc2iphi;
    
    
    bc2ieta = electronbc2ieta[j];
    bc2iphi = electronbc2iphi[j];
    bc2etacry = electronbc2etacry[j];
    bc2phicry = electronbc2phicry[j];
    fVals[64] = hasbc2 ? bc2ieta : 0.;
    fVals[65] = hasbc2 ? bc2iphi : 0.;
    fVals[66] = hasbc2 ? bc2ieta%5 : 0.;
    fVals[67] = hasbc2 ? bc2iphi%2 : 0.;
    fVals[68] = hasbc2 ? (TMath::Abs(bc2ieta)<=25)*(bc2ieta%25) + (TMath::Abs(bc2ieta)>25)*((bc2ieta-25*TMath::Abs(bc2ieta)/bc2ieta)%20) : 0.;
    fVals[69] = hasbc2 ? bc2iphi%20 : 0.;
    fVals[70] = hasbc2 ? bc2etacry : 0.;
    fVals[71] = hasbc2 ? bc2phicry : 0.;
        
    //Nb. of vertex
    fVals[72] = nVertex; 
    
    
  //   if( fabs( fVals[0]- 31.24) < 0.01){

//       for(int n=0; n<73; n++){
// 	cout<<"fVals["<<n<<"] " << fVals[n]<<endl;
//       }
//     }
    
    
  }else{
    
    fVals[0]  = electronscrawEnergy[j];
    fVals[1]  = electrone3x3[j]/electronscrawEnergy[j];
    fVals[2]  = electronsceta[j];
    fVals[3]  = electronscphi[j];
    fVals[4]  = electrone5x5[j]/electronscrawEnergy[j];
    fVals[5] = electronscetaWidth[j];
    fVals[6] = electronscphiWidth[j];
    //Nb. of vertex
    fVals[7] = nVertex; 
    
  }
  
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = electronscrawEnergy[j];
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = electronscrawEnergy[j] + electronscpreshowerEnergy[j]; 
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  
  //   if( isbarrel){
  //     return std::pair<double,double>(electronscenergy[j],0);
  //   }
  
  //cout<<"getting greder.." <<endl; 
  
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  
  //   if( isbarrel){
  //     cout<<"checkele " << greader->GetResponse(fVals) <<endl; 
  //     for(int j=0;j<15; j++){
  //       cout<<" fVals" << fVals[j]<<endl; 
  //     }
  //   }
    
  ///checked on 10000 events, no problem 
  ///cout<<"gett  "<< isbarrel<<" "<< electronscenergy[j]<<" "<< ecor <<" "<< ecorerr <<" ntpl " <<  electronenergyRegCorrected[j]<<" "<< electronenergyRegCorrectedError[j]<<endl; 
 //  if( fabs(ecor - electronenergyRegCorrected[j]) > 0.00003 || fabs(ecorerr - electronenergyRegCorrectedError[j]) > 0.00003 ){
//     cout<<"gett  "<< isbarrel<<" "<< evtNumber <<" "<< electronscenergy[j]<<" "<< ecor <<" "<< ecorerr <<" ntpl " <<  electronenergyRegCorrected[j]<<" "<< electronenergyRegCorrectedError[j]<< "diff " << ecor - electronenergyRegCorrected[j] <<" "<< ecorerr - electronenergyRegCorrectedError[j] <<" reldiff " << fabs(ecor - electronenergyRegCorrected[j])/ (0.5*(ecor + electronenergyRegCorrected[j])) <<endl; 
//   }
  
  return std::pair<double,double>(ecor,ecorerr);
  
  
}




std::pair<double,double> photonEnergyCorrector_CorrectedEnergyWithErrorv2(int j){
    
  
//   ///no correction for electron of tracker-deriven seed only ( not used in analysis anyway)
//   if( !electronecalDrivenSeed[j]){
//     return std::pair<double,double>(electronscenergy[j],0);
//   }
  

  getGapCoordinates(photonsceta[j],photonscphi[j]);
  
 

  Bool_t isbarrel = fabs(photonsceta[j]) < 1.48; 
  
  if( isbarrel){
        
    
    fVals[0]  = photonscrawEnergy[j];
    fVals[1]  = photone3x3[j]/photonscrawEnergy[j];
    fVals[2]  = photonsceta[j];
    fVals[3]  = photonscphi[j];
    fVals[4]  = photone5x5[j]/photonscrawEnergy[j];
    fVals[5] = photonhadronicOverEm[j];
    fVals[6] = photonscetaWidth[j];
    fVals[7] = photonscphiWidth[j];
        
    ///bc
    float bemax = photonbcseedeMax[j];
    float be2nd = photonbcseede2nd[j];
    float betop = photonbcseedeTop[j];
    float bebottom = photonbcseedeBottom[j];
    float beleft = photonbcseedeLeft[j];
    float beright = photonbcseedeRight[j];
    

    fVals[8]  = photonbcseedeta[j] - photonsceta[j]; 
    fVals[9]  = DeltaPhi(photonbcseedphi[j],photonscphi[j]);
    fVals[10]  = photonbcseedenergy[j]/ photonscrawEnergy[j];
    fVals[11]  = photonbcseede3x3[j] / photonbcseedenergy[j];
    fVals[12]  = photonbcseede5x5[j] / photonbcseedenergy[j];
    fVals[13] = sqrt(photonbcseedCovIEtaIEta[j]);
    fVals[14] = sqrt(photonbcseedCovIPhiIPhi[j]);
    fVals[15] = photonbcseedCovIEtaIPhi[j];
    fVals[16] = bemax/photonbcseedenergy[j];
    fVals[17] = log(be2nd/bemax);
    fVals[18] = log(betop/bemax);
    fVals[19] = log(bebottom/bemax);
    fVals[20] = log(beleft/bemax);
    fVals[21] = log(beright/bemax);
    fVals[22] = (betop-bebottom)/(betop+bebottom);
    fVals[23] = (beleft-beright)/(beleft+beright);

    

    ///bc2
    bool hasbc2 = photonbc2eMax[j] > 0; 
    float bc2emax = photonbc2eMax[j];
    float bc2e2nd = photonbc2e2nd[j];
    float bc2etop = photonbc2eTop[j];
    float bc2ebottom = photonbc2eBottom[j];
    float bc2eleft = photonbc2eLeft[j];
    float bc2eright = photonbc2eRight[j];
    
    fVals[24] = hasbc2 ? photonbc2eta[j] - photonsceta[j]: 0; 
    fVals[25] = hasbc2 ? DeltaPhi(photonbc2phi[j],photonscphi[j]): 0; 
    fVals[26] = hasbc2 ? photonbc2energy[j]/ photonscrawEnergy[j]: 0; 
    fVals[27] = hasbc2 ? photonbc2e3x3[j]/photonbc2energy[j]: 0; 
    fVals[28] = hasbc2 ? photonbc2e5x5[j]/photonbc2energy[j]: 0;
    fVals[29] = hasbc2 ? sqrt(photonbc2CovIEtaIEta[j]): 0;
    fVals[30] = hasbc2 ? sqrt(photonbc2CovIPhiIPhi[j]): 0; 
    
    
    fVals[31] = hasbc2 ? photonbc2CovIEtaIPhi[j]: 0; 
    ///bug ? 
    //fVals[31] = hasbc2 ? photonbcseedCovIEtaIPhi[j]: 0; 
    
    fVals[32] = hasbc2 ? bc2emax/photonbc2energy[j] : 0.;
    fVals[33] = hasbc2 ? log(bc2e2nd/bc2emax) : 0.;
    fVals[34] = hasbc2 ? log(bc2etop/bc2emax) : 0.;
    fVals[35] = hasbc2 ? log(bc2ebottom/bc2emax) : 0.;
    fVals[36] = hasbc2 ? log(bc2eleft/bc2emax) : 0.;
    fVals[37] = hasbc2 ? log(bc2eright/bc2emax) : 0.;
    fVals[38] = hasbc2 ? (bc2etop-bc2ebottom)/(bc2etop+bc2ebottom) : 0.;
    fVals[39] = hasbc2 ? (bc2eleft-bc2eright)/(bc2eleft+bc2eright) : 0.;
    
    //bclast
    bool hasbclast = photonbclastenergy[j]> 0; 
    fVals[40] = hasbclast ? (photonbclasteta[j]-photonsceta[j]) : 0.;
    fVals[41] = hasbclast ? DeltaPhi(photonbclastphi[j],photonscphi[j]) : 0.;
    fVals[42] = hasbclast ? photonbclastenergy[j] / photonscrawEnergy[j] : 0.;  
    fVals[43] = hasbclast ? photonbclaste3x3[j] / photonbclastenergy[j] :0.; 
    fVals[44] = hasbclast ? photonbclaste5x5[j] / photonbclastenergy[j] :0.; 
    fVals[45] = hasbclast ? sqrt(photonbclastCovIEtaIEta[j]): 0;
    fVals[46] = hasbclast ? sqrt(photonbclastCovIPhiIPhi[j]): 0; 
    fVals[47] = hasbclast ? photonbclastCovIEtaIPhi[j]: 0; 
    
    ///bclast2
    bool hasbclast2 = photonbclast2energy[j]> 0; 
    fVals[48] = hasbclast2 ? (photonbclast2eta[j]-photonsceta[j]) : 0.;
    fVals[49] = hasbclast2 ? DeltaPhi(photonbclast2phi[j],photonscphi[j]) : 0.;
    fVals[50] = hasbclast2 ? photonbclast2energy[j] / photonscrawEnergy[j] : 0.;  
    fVals[51] = hasbclast2 ? photonbclast2e3x3[j] / photonbclast2energy[j] :0.; 
    fVals[52] = hasbclast2 ? photonbclast2e5x5[j] / photonbclast2energy[j] :0.; 
    fVals[53] = hasbclast2 ? sqrt(photonbclast2CovIEtaIEta[j]): 0;
    fVals[54] = hasbclast2 ? sqrt(photonbclast2CovIPhiIPhi[j]): 0; 
    fVals[55] = hasbclast2 ? photonbclast2CovIEtaIPhi[j]: 0; 
    
    //local coordinates and crystal indices                                                                                                                                 
    //seed cluster                                                                                                                                                          
    float betacry, bphicry, bthetatilt, bphitilt;
    int bieta, biphi;
    
    bieta = photonbcseedieta[j];
    biphi = photonbcseediphi[j];
    betacry = photonbcseedetacry[j];
    bphicry = photonbcseedphicry[j];
    
    fVals[56] = bieta; //crystal ieta                                                                                                                                       
    fVals[57] = biphi; //crystal iphi                                                                                                                                       
    fVals[58] = bieta%5; //submodule boundary eta symmetry                                                                                                                  
    fVals[59] = biphi%2; //submodule boundary phi symmetry    
    fVals[60] = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry            
    fVals[61] = biphi%20; //module boundary phi symmetry                                                                                                                    
    fVals[62] = betacry; //local coordinates with respect to closest crystal center at nominal shower depth                                                                 
    fVals[63] = bphicry;

    //2nd cluster (meaningful gap corrections for converted photons)                                                                                                        
    float bc2etacry, bc2phicry, bc2thetatilt, bc2phitilt;
    int bc2ieta, bc2iphi;
    
    
    bc2ieta = photonbc2ieta[j];
    bc2iphi = photonbc2iphi[j];
    bc2etacry = photonbc2etacry[j];
    bc2phicry = photonbc2phicry[j];
    fVals[64] = hasbc2 ? bc2ieta : 0.;
    fVals[65] = hasbc2 ? bc2iphi : 0.;
    fVals[66] = hasbc2 ? bc2ieta%5 : 0.;
    fVals[67] = hasbc2 ? bc2iphi%2 : 0.;
    fVals[68] = hasbc2 ? (TMath::Abs(bc2ieta)<=25)*(bc2ieta%25) + (TMath::Abs(bc2ieta)>25)*((bc2ieta-25*TMath::Abs(bc2ieta)/bc2ieta)%20) : 0.;
    fVals[69] = hasbc2 ? bc2iphi%20 : 0.;
    fVals[70] = hasbc2 ? bc2etacry : 0.;
    fVals[71] = hasbc2 ? bc2phicry : 0.;
        
    //Nb. of vertex
    fVals[72] = nVertex; 
    
    
    
  }else{
    
    fVals[0]  = photonscrawEnergy[j];
    fVals[1]  = photone3x3[j]/photonscrawEnergy[j];
    fVals[2]  = photonsceta[j];
    fVals[3]  = photonscphi[j];
    fVals[4]  = photone5x5[j]/photonscrawEnergy[j];
    fVals[5] = photonscetaWidth[j];
    fVals[6] = photonscphiWidth[j];
    //Nb. of vertex
    fVals[7] = nVertex; 
    
//     if( evtNumber == 369579111 || evtNumber == 372222387 ){
//       for(int n=0; n<8; n++){
// 	cout<<"phteefVals["<<n<<"] " << fVals[n]<<endl;
//       }
//     }
    
    
  }
  
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = photonscrawEnergy[j];
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = photonscrawEnergy[j] + photonscpreshowerEnergy[j]; 
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  
  //   if( isbarrel){
  //     return std::pair<double,double>(photonscenergy[j],0);
  //   }
  
  //cout<<"getting greder.." <<endl; 
  
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  
  //   if( isbarrel){
  //     cout<<"checkele " << greader->GetResponse(fVals) <<endl; 
  //     for(int j=0;j<15; j++){
  //       cout<<" fVals" << fVals[j]<<endl; 
  //     }
  //   }
    
  ///checked on 10000 events, no problem 
  //cout<<"gett  "<< isbarrel<<" "<< photonscenergy[j]<<" "<< ecor <<" "<< ecorerr <<" ntpl " <<  photonenergyRegCorrected[j]<<" "<< photonenergyRegCorrectedError[j]<<endl; 

//   if( fabs(ecor - photonenergyRegCorrected[j]) > 0.00003 || fabs(ecorerr - photonenergyRegCorrectedError[j]) > 0.00003 ){
//     cout<<"gett  "<< isbarrel<<" "<< evtNumber <<" "<< photonscenergy[j]<<" "<< ecor <<" "<< ecorerr <<" ntpl " <<  photonenergyRegCorrected[j]<<" "<< photonenergyRegCorrectedError[j]<< "diff " << ecor - photonenergyRegCorrected[j] <<" "<< ecorerr - photonenergyRegCorrectedError[j] <<" reldiff " << fabs(ecor - photonenergyRegCorrected[j])/ (0.5*(ecor + photonenergyRegCorrected[j])) <<endl; 
//}
  
  
  return std::pair<double,double>(ecor,ecorerr);
  
  
}
