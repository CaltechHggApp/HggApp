#include "SusyHggSelectorSimple.hh"

#include "HggSelector.hh"

#define VERBOSE 0

SusyHggSelectorSimple::SusyHggSelectorSimple(std::vector<std::string> fNames, std::string treeName,std::string outputFile):
  SusyHggSelector(fNames,treeName,outputFile)
{
  min_ptgg=0;
}

void SusyHggSelectorSimple::processConfig(ReadConfig& cfg) {
}

void SusyHggSelectorSimple::setSmearingCFG(std::string smearCFG) {
  if(smearer!=0) delete smearer;
  smearer = new HggEnergyScale(smearCFG);
}

void SusyHggSelectorSimple::setScaleingCFG(std::string scaleCFG) {
  if(scaleer!=0) delete scaleer;
  scaleer = new HggEnergyScale(scaleCFG);
}

void SusyHggSelectorSimple::clear() {
  mgg=-1;
  ptgg=-1;
  
  pho1_pt=-1;
  pho1_eta=0;
  pho1_phi=0;
  pho1_r9=-1;
  
  pho2_pt=-1;
  pho2_eta=0;
  pho2_phi=0;
  pho2_r9=-1;

  hem1_pt=0;
  hem2_pt=0;

  mgg_def=-1;
  ptgg_def=-1;
  
  pho1_pt_def=-1;
  pho1_eta_def=0;
  pho1_phi_def=0;
  pho1_r9_def=-1;
  
  pho2_pt_def=-1;
  pho2_eta_def=0;
  pho2_phi_def=0;
  pho2_r9_def=-1;

  hem1_def_pt=0;
  hem2_def_pt=0;

  ele1_pt=0;
  mu1_pt=0;

  nJ=0;
  MR=0;
  Rsq=0;
  t1Rsq=0;

  MR_def=0;
  Rsq_def=0;
  t1Rsq_def=0;

  nSusyPart=0;
  for(int i=0;i<20;i++) {
    mSusyPart[i]=0;
    idSusyPart[i]=0;
  }

  m22=0;
  m23=0;
  m24=0;
  m25=0;

  for(int i=0;i<kMaxJets;i++) {
    indexJet[i]=-1;
    ptJet[i]=0;
    etaJet[i]=999;
    phiJet[i]=999;
    energyJet[i]=0;
    hemJet[i]=0;
  }
}


void SusyHggSelectorSimple::processEntry(Long64_t iEntry) {
  if(nVtx==0) return;
  
  TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);

  puWeight=pileupWeight;
  runNum  = runNumber;
  lumiSec = lumiBlock;
  evtNum  = evtNumber;

  outRho = rho;

    for(int i=0;i<nPho_;i++) {
      auto photon = &(Photons_->at(i));
      if(isMC) { 
	if(smearer) {//smear the photons
	  photon->scaledEnergy = photon->correctedEnergy;
	  photon->scaledEnergyError = photon->correctedEnergyError;
	  std::pair<float,float> dE = smearer->getDEoE(*photon,0);
	  photon->dEoE = dE.first;
	  photon->dEoEErr = dE.second;
	  
	  assert(photon!=0);
	  HggSelector::smearPhoton(photon,0);
	}else{ //no smearer
	  photon->finalEnergy = photon->correctedEnergy;
	  photon->finalEnergyError = photon->correctedEnergyError;
	}
	//std::cout <<">>>> "<< photon->correctedEnergy << "  " << photon->finalEnergy << "  " << photon->correctedEnergyError << "  " << photon->finalEnergyError  << "  " << photon->dEoE << std::endl;
      }else{ //data
	if(scaler) {
	  std::pair<float,float> dE = energyScale->getDEoE(*photon,runNumber);
	  photon->dEoE    = dE.first;
	  photon->scaledEnergy = photon->correctedEnergy*(photon->dEoE);
	  photon->scaledEnergyError = photon->correctedEnergyError*((photon->dEoE+photon->dEoEErr));
	  photon->finalEnergy = photon->scaledEnergy;
	  photon->finalEnergyError = photon->scaledEnergyError;
	}else{//no scaler
	  photon->finalEnergy = photon->correctedEnergy;
	  photon->finalEnergyError = photon->correctedEnergyError;
	}
      }
    }

  //identify the two highest pT photons
  const int NPho=2;
  float pho_sum_pt=0;
  float pho_sum_pt_def=0;
  float M_max=-1;

  int selected_photons[NPho];
  for(int i=0;i<NPho;i++) {
    selected_photons[i]=0;
  }
  int selected_photons_def[NPho];
  for(int i=0;i<NPho;i++) {
    selected_photons_def[i]=0;
  }
  for(int iPho=0; iPho < nPho_;iPho++) {
    auto photon1 =&(Photons_->at(iPho));
    //std::cout << photon1->correctedEnergy << "  " << photon1->finalEnergy << "  " << photon1->correctedEnergyError << "  " << photon1->finalEnergyError  << "  " << photon1->dEoE << std::endl;
    //veto photons outside acceptance
    if(VERBOSE) std::cout << photon1->SC.eta << std::endl;
    if(fabs(photon1->SC.eta) > max_pho_eta) continue;
    if(fabs(photon1->SC.eta) > 1.4442 && fabs(photon1->SC.eta) < 1.566) continue; //gap photons

    float pt1 = photon1->energy/cosh(photon1->eta);
    //float pt1 = photon1->p4FromVtx(vtx,photon1->finalEnergy).Pt();
    if(!optimize) {
      if(pt1 < min_pho2_pt) continue; //not hard enough
      if(VERBOSE) std::cout << "\tpass pt cut" << std::endl;
      if( !photonID.passID(*photon1,StandardPhotonID::kLoose,rho) ) continue; //no photon ID
      if( !photonID.passIso(*photon1,StandardPhotonID::kLoose,rho) ) continue; //no photon Iso

      if(VERBOSE) std::cout << "\tpass ID cut" << std::endl;

      std::bitset<5> id_res = photonID.cutResults(*photon1,StandardPhotonID::kLoose,rho);
      //if(id_res[3] || id_res[2]) continue; //fails the charged hadron or neutral hadron ID
      
      if(!isDY && photonMatchedElectron[iPho]) continue; //conversion-safe electron veto
      if(VERBOSE) std::cout << "\tpass CSEV cut" << std::endl;
    }
    
    TLorentzVector p4_pho1 = photon1->p4FromVtx(vtx,photon1->finalEnergy);
    TLorentzVector p4_pho1_def = photon1->p4FromVtx(vtx,photon1->energy);

    for(int jPho=iPho+1;jPho<nPho_;jPho++) {//second photon
      auto photon2 = &(Photons_->at(jPho));
      //veto photons outside acceptance
      if(fabs(photon2->SC.eta) > max_pho_eta) continue;
      if(fabs(photon2->SC.eta) > 1.4442 && fabs(photon1->SC.eta) < 1.566) continue; //gap photons
      
      float pt2 = photon2->energy/cosh(photon2->eta);
      //float pt2 = photon2->p4FromVtx(vtx,photon2->finalEnergy).Pt();
      if(!optimize) {
	if(pt2 < min_pho2_pt) continue; //not hard enough
	if(pt1 <min_pho1_pt && pt2 < min_pho1_pt) continue; //neither is hard enough
	if(VERBOSE) std::cout << "\tpass 40/25 cut" << std::endl;

	if( !photonID.passID(*photon2,StandardPhotonID::kLoose,rho) ) continue; //no photon ID	
	if( !photonID.passIso(*photon2,StandardPhotonID::kLoose,rho) ) continue; //no photon Iso

	std::bitset<5> id_res = photonID.cutResults(*photon2,StandardPhotonID::kLoose,rho);
	//if(id_res[3] || id_res[2]) continue; //fails the charged hadron or neutral hadron ID
	
	if(photonMatchedElectron[jPho]) continue; //conversion-safe electron veto
      }
      TLorentzVector p4_pho2 = photon2->p4FromVtx(vtx,photon2->finalEnergy);
      TLorentzVector p4_pho2_def = photon2->p4FromVtx(vtx,photon2->energy);
      float M = (p4_pho1+p4_pho2).M();
      if(VERBOSE) std::cout << "\tM: " << M  << std::endl;
      if( (p4_pho1.Pt() + p4_pho2.Pt()) > pho_sum_pt) {
	selected_photons[0]=(pt1>pt2 ? iPho : jPho);
	selected_photons[1]=(pt1>pt2 ? jPho : iPho);
	pho_sum_pt = (p4_pho1.Pt() + p4_pho2.Pt());
	if(VERBOSE) std::cout << "\tsumPT: " << pho_sum_pt << std::endl;
      }
      float M_def = (p4_pho1_def+p4_pho2_def).M();
      if(VERBOSE) std::cout << "\tM: " << M  << std::endl;
      if( (p4_pho1_def.Pt() + p4_pho2_def.Pt()) > pho_sum_pt_def) {
	selected_photons_def[0]=(pt1>pt2 ? iPho : jPho);
	selected_photons_def[1]=(pt1>pt2 ? jPho : iPho);
	pho_sum_pt_def = (p4_pho1_def.Pt() + p4_pho2_def.Pt());
	if(VERBOSE) std::cout << "\tsumPT: " << pho_sum_pt << std::endl;
      }
    }
  }
  if(selected_photons[0]==selected_photons[1] && selected_photons_def[0]==selected_photons_def[1]) return; //no selected photons

  //std::cout << "pho pts: " << pho_pts[0] << "  "<< pho_pts[1] << std::endl;
  //pho_pts is now the sorted pts of the highest pt photons in the event

  //selected_photons[] now contains the (sorted) reference to the two highest pT photons

  TLorentzVector pho1_p4 = Photons_->at(selected_photons[0]).p4FromVtx(vtx,Photons_->at(selected_photons[0]).finalEnergy);
  TLorentzVector pho2_p4 = Photons_->at(selected_photons[1]).p4FromVtx(vtx,Photons_->at(selected_photons[1]).finalEnergy);

  TLorentzVector pho1_p4_def = Photons_->at(selected_photons_def[0]).p4FromVtx(vtx,Photons_->at(selected_photons_def[0]).energy);
  TLorentzVector pho2_p4_def = Photons_->at(selected_photons_def[1]).p4FromVtx(vtx,Photons_->at(selected_photons_def[1]).energy);

  TLorentzVector gg_p4 = pho1_p4+pho2_p4;
  TLorentzVector gg_p4_def = pho1_p4_def+pho2_p4_def;

  mgg = gg_p4.M();
  mgg_def = gg_p4_def.M();
  //std::cout << "mgg: " << mgg << std::endl;
  ptgg = gg_p4.Pt();
  ptgg_def = gg_p4_def.Pt();

  if(VERBOSE) std::cout << ">>>>>>>>>>>>>" << selected_photons[0] << " " << selected_photons[1] << " " << mgg << "  " << ptgg << std::endl;
  if(!optimize) {
    if(mgg < min_mgg && mgg_def < min_mgg) return; //mgg too low
    if(ptgg < min_ptgg && ptgg_def < min_ptgg) return;
  }
  etagg = gg_p4.Eta();
  phigg = gg_p4.Phi();
  pho1_pt = pho1_p4.Pt();
  pho1_eta = pho1_p4.Eta();
  pho1_sc_eta = Photons_->at(selected_photons[0]).SC.eta;
  pho1_phi = pho1_p4.Phi();
  pho1_r9 = Photons_->at(selected_photons[0]).SC.r9;
  //if(isMC)pho1_seoe = Photons_->at(selected_photons[0]).correctedEnergyError/Photons_->at(selected_photons[0]).correctedEnergy;
  pho1_seoe = Photons_->at(selected_photons[0]).finalEnergyError/Photons_->at(selected_photons[0]).finalEnergy;

  std::bitset<5> id_res_pho1 = photonID.cutResults(Photons_->at(selected_photons[0]),StandardPhotonID::kLoose,rho);

  pho1_pass_id  = !(id_res_pho1[0] || id_res_pho1[1]);
  pho1_pass_iso = !(id_res_pho1[2] || id_res_pho1[3] || id_res_pho1[4]);


  pho1_sieie   = Photons_->at(selected_photons[0]).SC.sigmaIEtaIEta;
  pho1_HE      = Photons_->at(selected_photons[0]).HoverE;
  pho1_charged = Photons_->at(selected_photons[0]).dr03ChargedHadronPFIso[0];
  pho1_neutral = Photons_->at(selected_photons[0]).dr03NeutralHadronPFIso;
  pho1_photon  = Photons_->at(selected_photons[0]).dr03PhotonPFIso;
  pho1_eleveto = photonMatchedElectron[selected_photons[0]];
  pho1_genMatch = (Photons_->at(selected_photons[0]).genMatch.index!=-1);

  pho1_energyGen = (pho1_genMatch ? Photons_->at(selected_photons[0]).genMatch.energy : -1);

  etagg_def = gg_p4_def.Eta();
  phigg_def = gg_p4_def.Phi();
  pho1_def_pt = pho1_p4_def.Pt();
  pho1_def_eta = pho1_p4_def.Eta();
  pho1_def_sc_eta = Photons_->at(selected_photons_def[0]).SC.eta;
  pho1_def_phi = pho1_p4_def.Phi();
  pho1_def_r9 = Photons_->at(selected_photons_def[0]).SC.r9;
  //if(isMC)pho1_def_seoe = Photons_->at(selected_photons[0]).correctedEnergyError/Photons_->at(selected_photons[0]).correctedEnergy;
  pho1_def_seoe = 0;

  std::bitset<5> id_res_pho1_def = photonID.cutResults(Photons_->at(selected_photons_def[0]),StandardPhotonID::kLoose,rho);

  pho1_def_pass_id  = !(id_res_pho1_def[0] || id_res_pho1_def[1]);
  pho1_def_pass_iso = !(id_res_pho1_def[2] || id_res_pho1_def[3] || id_res_pho1_def[4]);


  pho1_def_sieie = Photons_->at(selected_photons_def[0]).SC.sigmaIEtaIEta;
  pho1_def_HE = Photons_->at(selected_photons_def[0]).HoverE;
  pho1_def_charged = Photons_->at(selected_photons_def[0]).dr03ChargedHadronPFIso[0];
  pho1_def_neutral = Photons_->at(selected_photons_def[0]).dr03NeutralHadronPFIso;
  pho1_def_photon = Photons_->at(selected_photons_def[0]).dr03PhotonPFIso;
  pho1_def_eleveto = photonMatchedElectron[selected_photons_def[0]];
  pho1_def_genMatch = (Photons_->at(selected_photons_def[0]).genMatch.index!=-1);

  pho1_def_energyGen = (pho1_def_genMatch ? Photons_->at(selected_photons_def[0]).genMatch.energy : -1);

  pho2_pt = pho2_p4.Pt();
  pho2_eta = pho2_p4.Eta();
  pho2_sc_eta = Photons_->at(selected_photons[1]).SC.eta;
  pho2_phi = pho2_p4.Phi();
  pho2_r9 = Photons_->at(selected_photons[1]).SC.r9;
  //if(isMC) pho2_seoe = Photons_->at(selected_photons[1]).correctedEnergyError/Photons_->at(selected_photons[1]).correctedEnergy;
  pho2_seoe = Photons_->at(selected_photons[1]).finalEnergyError/Photons_->at(selected_photons[1]).finalEnergy;

  std::bitset<5> id_res_pho2 = photonID.cutResults(Photons_->at(selected_photons[1]),StandardPhotonID::kLoose,rho);

  pho2_pass_id  = !(id_res_pho2[0] || id_res_pho2[1]);
  pho2_pass_iso = !(id_res_pho2[2] || id_res_pho2[3] || id_res_pho2[4]);

  pho2_sieie   = Photons_->at(selected_photons[1]).SC.sigmaIEtaIEta;
  pho2_HE      = Photons_->at(selected_photons[1]).HoverE;
  pho2_charged = Photons_->at(selected_photons[1]).dr03ChargedHadronPFIso[0];
  pho2_neutral = Photons_->at(selected_photons[1]).dr03NeutralHadronPFIso;
  pho2_photon  = Photons_->at(selected_photons[1]).dr03PhotonPFIso;
  pho2_eleveto = photonMatchedElectron[selected_photons[1]];
  pho2_genMatch = (Photons_->at(selected_photons[1]).genMatch.index!=-1);

  pho2_energyGen = (pho2_genMatch ? Photons_->at(selected_photons[1]).genMatch.energy : -1);


  pho2_def_pt = pho2_p4_def.Pt();
  pho2_def_eta = pho2_p4_def.Eta();
  pho2_def_sc_eta = Photons_->at(selected_photons_def[1]).SC.eta;
  pho2_def_phi = pho2_p4_def.Phi();
  pho2_def_r9 = Photons_->at(selected_photons_def[1]).SC.r9;
  //if(isMC) pho2_def_seoe = Photons_->at(selected_photons_def[1]).correctedEnergyError/Photons_->at(selected_photons_def[1]).correctedEnergy;
  pho2_def_seoe = 0; 

  std::bitset<5> id_res_pho2_def = photonID.cutResults(Photons_->at(selected_photons_def[1]),StandardPhotonID::kLoose,rho);

  pho2_def_pass_id  = !(id_res_pho2_def[0] || id_res_pho2_def[1]);
  pho2_def_pass_iso = !(id_res_pho2_def[2] || id_res_pho2_def[3] || id_res_pho2_def[4]);

  pho2_def_sieie = Photons_->at(selected_photons_def[1]).SC.sigmaIEtaIEta;
  pho2_def_HE = Photons_->at(selected_photons_def[1]).HoverE;
  pho2_def_charged = Photons_->at(selected_photons_def[1]).dr03ChargedHadronPFIso[0];
  pho2_def_neutral = Photons_->at(selected_photons_def[1]).dr03NeutralHadronPFIso;
  pho2_def_photon = Photons_->at(selected_photons_def[1]).dr03PhotonPFIso;
  pho2_def_eleveto = photonMatchedElectron[selected_photons_def[1]];
  pho2_def_genMatch = (Photons_->at(selected_photons_def[1]).genMatch.index!=-1);

  pho2_def_energyGen = (pho2_def_genMatch ? Photons_->at(selected_photons_def[1]).genMatch.energy : -1);

  //clear the highest_csv info
  highest_csv = -1000;
  highest_csv_pt=0;
  highest_csv_eta=0;
  highest_csv_phi=0;

  highest_csv_up = -1000;
  highest_csv_pt_up=0;
  highest_csv_eta_up=0;
  highest_csv_phi_up=0;

  highest_csv_down = -1000;
  highest_csv_pt_down=0;
  highest_csv_eta_down=0;
  highest_csv_phi_down=0;

  //clear the highest_csv info
  second_csv = -1000;
  second_csv_pt=0;
  second_csv_eta=0;
  second_csv_phi=0;

  second_csv_up = -1000;
  second_csv_pt_up=0;
  second_csv_eta_up=0;
  second_csv_phi_up=0;

  second_csv_down = -1000;
  second_csv_pt_down=0;
  second_csv_eta_down=0;
  second_csv_phi_down=0;

  //loop over the jets
  std::vector<TLorentzVector> selectedJets;
  std::vector<TLorentzVector> selectedJets_up;
  std::vector<TLorentzVector> selectedJets_down;

  selectJets(&selectedJets,0);
  selectJets(&selectedJets_up,1);
  selectJets(&selectedJets_down,-1);

  auto calcVars = [](float& HT,float& MHT,float& MHTphi,std::vector<TLorentzVector>& jets,TLorentzVector& pho1, TLorentzVector& pho2) -> void
    {
      HT=0;
      TLorentzVector MHTvec; MHTvec.SetPtEtaPhiM(0.,0.,0.,0.);
      for(auto jet : jets) {
	HT+=jet.Pt();
	MHTvec-=jet;
      }
      HT+=pho1.Pt();
      HT+=pho2.Pt();
      MHTvec-=pho1;
      MHTvec-=pho2;
      
      MHT = MHTvec.Pt();
      MHTphi = MHTvec.Phi();
    };

  calcVars(HT,MHT,MHTphi,selectedJets,pho1_p4,pho2_p4);
  calcVars(HT_up,MHT_up,MHTphi_up,selectedJets_up,pho1_p4,pho2_p4);
  calcVars(HT_down,MHT_down,MHTphi_down,selectedJets_down,pho1_p4,pho2_p4);

  TVector3 met;
  met.SetPtEtaPhi(pfMet,0.,pfMetPhi);
  TVector3 t1met;
  t1met.SetPtEtaPhi(type1PfMet,0.,type1PfMetPhi);
  
  MET = pfMet;
  METphi = pfMetPhi;

  t1MET = type1PfMet;
  t1METphi = type1PfMetPhi;

  nJ = selectedJets.size();
  nJ_up = selectedJets_up.size();
  nJ_down = selectedJets_down.size();

  if(gg_p4.M()>0) { //there is a regression pair
    selectedJets.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
    selectedJets_up.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
    selectedJets_down.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
    if(!optimize) {
      if(selectedJets.size()<2) {
	if(VERBOSE) std::cout << "NOT ENOUGH JETS" << std::endl;
	return;
      }
    }

    TLorentzVector top_b;
    TLorentzVector second_b;

    if(selectedJets.size()>=2) {
      TLorentzVector h1(0,0,0,0),h2(0,0,0,0);
      std::vector<int> hemAssign;
      try{
	RazorVariables::CombineJets(selectedJets,h1,h2,&hemAssign);
      }catch(TooManyJets& e) {
	std::cout << "TOO MANY JETS IN EVENT" << std::endl;
	return;
      }catch(TooFewJets& e) {
	std::cout << "shouldn't happen......" << std::endl;
	throw e;
      }

      for(int i=0;i<selectedJets.size();i++) {
	//fill the jet info
	indexJet[i] = i;
	ptJet[i] = selectedJets[i].Pt();
	etaJet[i] = selectedJets[i].Eta();
	phiJet[i] = selectedJets[i].Phi();
	energyJet[i] = selectedJets[i].E();      

	try{
	  corrUpJet[i] = jecReader.getCorrection(ptJet[i],etaJet[i],JECUReader::kUp);
	  corrDownJet[i] = jecReader.getCorrection(ptJet[i],etaJet[i],JECUReader::kDown);
	}catch(std::runtime_error& e){
	  corrUpJet[i]=0;
	  corrDownJet[i]=0;
	}
	hemJet[i] = hemAssign.at(i);
      }

      hemgg = hemAssign.at(hemAssign.size()-1);
      if(h1.Pt()==0) {
	if(VERBOSE) std::cout << "HEMISPHERE 1 had 0 pT!" << std::endl;
	return;
      }
      hem1_pt  = h1.Pt();
      hem1_eta = h1.Eta();
      hem1_phi = h1.Phi();
      hem1_M   = h1.M();
      
      hem2_pt  = h2.Pt();
      hem2_eta = h2.Eta();
      hem2_phi = h2.Phi();
      hem2_M   = h2.M();
      
      MR = RazorVariables::CalcGammaMRstar(h1,h2);

      double MTR = RazorVariables::CalcMTR(h1,h2,met);
      double t1MTR = RazorVariables::CalcMTR(h1,h2,t1met);

      Rsq = MTR/MR;
      Rsq=Rsq*Rsq; //square it!

      t1Rsq = t1MTR/MR;
      t1Rsq=t1Rsq*t1Rsq; //square it!
      
    
    }
    
    if(selectedJets_up.size()>=2) {
      TLorentzVector h1(0,0,0,0),h2(0,0,0,0);
      try{
	RazorVariables::CombineJets(selectedJets_up,h1,h2);
      }catch(TooManyJets& e) {
	std::cout << "TOO MANY JETS IN EVENT" << std::endl;
	return;
      }catch(TooFewJets& e) {
	std::cout << "shouldn't happen......" << std::endl;
	throw e;
      }
      
      if(h1.Pt()!=0) {
	
	hem1_pt_up  = h1.Pt();
	hem1_eta_up = h1.Eta();
	hem1_phi_up = h1.Phi();
	hem1_M_up   = h1.M();
	
	hem2_pt_up  = h2.Pt();
	hem2_eta_up = h2.Eta();
	hem2_phi_up = h2.Phi();
	hem2_M_up   = h2.M();
	
	MR_up = RazorVariables::CalcGammaMRstar(h1,h2);
	
	double MTR = RazorVariables::CalcMTR(h1,h2,met);
	double t1MTR = RazorVariables::CalcMTR(h1,h2,t1met);
	
	Rsq_up = MTR/MR_up;
	Rsq_up=Rsq_up*Rsq_up; //square it!
	
	t1Rsq_up = t1MTR/MR_up;
	t1Rsq_up=t1Rsq_up*t1Rsq_up; //square it!
      }else{
	t1Rsq_up=-1;
	Rsq_up=-1;
	MR_up=-1;
      }
    }
    
    if(selectedJets_down.size()>=2) {
      TLorentzVector h1(0,0,0,0),h2(0,0,0,0);
      try{
	RazorVariables::CombineJets(selectedJets_down,h1,h2);
      }catch(TooManyJets& e) {
	std::cout << "TOO MANY JETS IN EVENT" << std::endl;
	return;
      }catch(TooFewJets& e) {
       std::cout << "shouldn't happen......" << std::endl;
       throw e;
      }

      if(h1.Pt()!=0) {

	hem1_pt_down  = h1.Pt();
	hem1_eta_down = h1.Eta();
	hem1_phi_down = h1.Phi();
	hem1_M_down   = h1.M();
	
	hem2_pt_down  = h2.Pt();
	hem2_eta_down = h2.Eta();
	hem2_phi_down = h2.Phi();
	hem2_M_down   = h2.M();
      
	MR_down = RazorVariables::CalcGammaMRstar(h1,h2);
	
	double MTR = RazorVariables::CalcMTR(h1,h2,met);
	double t1MTR = RazorVariables::CalcMTR(h1,h2,t1met);
	
	t1Rsq_down = t1MTR/MR_down;
	t1Rsq_down=t1Rsq_down*t1Rsq_down; //square it!
      }else{
	t1Rsq_down=-1;
	Rsq_down=-1;
	MR_down=-1;
      }
    }

    //remove the regression Higgs
    selectedJets.pop_back();
    selectedJets_up.pop_back();
    selectedJets_down.pop_back();
    
  }else{
    MR=0;
    Rsq=0;
    t1Rsq=0;
    MR_up=0;
    Rsq_up=0;
    t1Rsq_up=0;
    MR_down=0;
    Rsq_down=0;
    t1Rsq_down=0;
  }


  if(gg_p4_def.M()>0) { //there is a non-regression pair
    selectedJets.push_back(gg_p4_def); //ensure the 'Higgs' vector is in there explicitly
    selectedJets_up.push_back(gg_p4_def); //ensure the 'Higgs' vector is in there explicitly
    selectedJets_down.push_back(gg_p4_def); //ensure the 'Higgs' vector is in there explicitly
    if(!optimize) {
      if(selectedJets.size()<2) {
	if(VERBOSE) std::cout << "NOT ENOUGH JETS" << std::endl;
	return;
      }
    }

    TLorentzVector top_b;
    TLorentzVector second_b;

    if(selectedJets.size()>=2) {
      TLorentzVector h1(0,0,0,0),h2(0,0,0,0);
      std::vector<int> hemAssign;
      try{
	RazorVariables::CombineJets(selectedJets,h1,h2,&hemAssign);
      }catch(TooManyJets& e) {
	std::cout << "TOO MANY JETS IN EVENT" << std::endl;
	return;
      }catch(TooFewJets& e) {
	std::cout << "shouldn't happen......" << std::endl;
	throw e;
      }

      for(int i=0;i<selectedJets.size();i++) {
	//fill the jet info
	indexJet[i] = i;
	ptJet[i] = selectedJets[i].Pt();
	etaJet[i] = selectedJets[i].Eta();
	phiJet[i] = selectedJets[i].Phi();
	energyJet[i] = selectedJets[i].E();      

	try{
	  corrUpJet[i] = jecReader.getCorrection(ptJet[i],etaJet[i],JECUReader::kUp);
	  corrDownJet[i] = jecReader.getCorrection(ptJet[i],etaJet[i],JECUReader::kDown);
	}catch(std::runtime_error& e){
	  corrUpJet[i]=0;
	  corrDownJet[i]=0;
	}
	hemJet_def[i] = hemAssign.at(i);
      }

      hemgg_def = hemAssign.at(hemAssign.size()-1);
      if(h1.Pt()>0) {

	hem1_def_pt  = h1.Pt();
	hem1_def_eta = h1.Eta();
	hem1_def_phi = h1.Phi();
	hem1_def_M   = h1.M();
      
	hem2_def_pt  = h2.Pt();
	hem2_def_eta = h2.Eta();
	hem2_def_phi = h2.Phi();
	hem2_def_M   = h2.M();
      
	MR_def = RazorVariables::CalcGammaMRstar(h1,h2);
	
	double MTR = RazorVariables::CalcMTR(h1,h2,met);
	double t1MTR = RazorVariables::CalcMTR(h1,h2,t1met);
	
	Rsq_def = MTR/MR_def;
	Rsq_def=Rsq_def*Rsq_def; //square it!
	
	t1Rsq_def = t1MTR/MR_def;
	t1Rsq_def=t1Rsq_def*t1Rsq_def; //square it!
      }else{
	MR_def=0;
	Rsq_def=0;
	t1Rsq_def=0;
      }
    
    }
    
    if(selectedJets_up.size()>=2) {
      TLorentzVector h1(0,0,0,0),h2(0,0,0,0);
      try{
	RazorVariables::CombineJets(selectedJets_up,h1,h2);
      }catch(TooManyJets& e) {
	std::cout << "TOO MANY JETS IN EVENT" << std::endl;
	return;
      }catch(TooFewJets& e) {
	std::cout << "shouldn't happen......" << std::endl;
	throw e;
      }
      
      if(h1.Pt()!=0) {
	
	hem1_def_pt_up  = h1.Pt();
	hem1_def_eta_up = h1.Eta();
	hem1_def_phi_up = h1.Phi();
	hem1_def_M_up   = h1.M();
	
	hem2_def_pt_up  = h2.Pt();
	hem2_def_eta_up = h2.Eta();
	hem2_def_phi_up = h2.Phi();
	hem2_def_M_up   = h2.M();
	
	MR_up_def = RazorVariables::CalcGammaMRstar(h1,h2);
	
	double MTR = RazorVariables::CalcMTR(h1,h2,met);
	double t1MTR = RazorVariables::CalcMTR(h1,h2,t1met);
	
	Rsq_up_def = MTR/MR_up_def;
	Rsq_up_def=Rsq_up_def*Rsq_up_def; //square it!
	
	t1Rsq_up_def = t1MTR/MR_up_def;
	t1Rsq_up_def=t1Rsq_up_def*t1Rsq_up_def; //square it!
      }else{
	t1Rsq_up_def=-1;
	Rsq_up_def=-1;
	MR_up_def=-1;
      }
    }
    
    if(selectedJets_down.size()>=2) {
      TLorentzVector h1(0,0,0,0),h2(0,0,0,0);
      try{
	RazorVariables::CombineJets(selectedJets_down,h1,h2);
      }catch(TooManyJets& e) {
	std::cout << "TOO MANY JETS IN EVENT" << std::endl;
	return;
      }catch(TooFewJets& e) {
       std::cout << "shouldn't happen......" << std::endl;
       throw e;
      }

      if(h1.Pt()!=0) {

	hem1_def_pt_down  = h1.Pt();
	hem1_def_eta_down = h1.Eta();
	hem1_def_phi_down = h1.Phi();
	hem1_def_M_down   = h1.M();
	
	hem2_def_pt_down  = h2.Pt();
	hem2_def_eta_down = h2.Eta();
	hem2_def_phi_down = h2.Phi();
	hem2_def_M_down   = h2.M();
      
	MR_down_def = RazorVariables::CalcGammaMRstar(h1,h2);
	
	double MTR = RazorVariables::CalcMTR(h1,h2,met);
	double t1MTR = RazorVariables::CalcMTR(h1,h2,t1met);
	
	t1Rsq_down_def = t1MTR/MR_down_def;
	t1Rsq_down_def=t1Rsq_down_def*t1Rsq_down_def; //square it!
      }else{
	t1Rsq_down_def=-1;
	Rsq_down_def=-1;
	MR_down_def=-1;
      }
    }
  }else{
    MR_def=0;
    Rsq_def=0;
    t1Rsq_def=0;
    MR_up_def=0;
    Rsq_up_def=0;
    t1Rsq_up_def=0;
    MR_down_def=0;
    Rsq_down_def=0;
    t1Rsq_down_def=0;
  }


  //find leptons
  int highest_ele=-1;
  float highest_ele_pt=-1;
  for(int iEle=0; iEle<nEle_; iEle++) {
    auto electron = Electrons_->at(iEle);
    if(! electronID.passAll(electron,StandardElectronID::kMedium) ) continue;
    
    if(electron.pt > highest_ele_pt) {
      highest_ele = iEle;
      highest_ele_pt = electron.pt;
    }
  }
  if(highest_ele!=-1) {
    ele1_pt = Electrons_->at(highest_ele).pt;
    ele1_eta = Electrons_->at(highest_ele).eta;
    ele1_phi = Electrons_->at(highest_ele).phi;
  }


  int highest_mu=-1;
  float highest_mu_pt=-1;
  for(int iMu=0; iMu<nMu_; iMu++) {
    auto mu = Muons_->at(iMu);
    
    if(! mu.isTightMuon ) continue;
    
    if(mu.pt > highest_mu_pt) {
      highest_mu = iMu;
      highest_mu_pt = mu.pt;
    }
  }

  if(highest_mu!=-1) {
    mu1_pt = Muons_->at(highest_mu).pt;
    mu1_eta = Muons_->at(highest_mu).eta;
    mu1_phi = Muons_->at(highest_mu).phi;
  }
  
  if(isMC) fillGenTruth();
  outTree->Fill(); //we need to fill ourselves

  

}



