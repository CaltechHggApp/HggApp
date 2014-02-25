#include "SusyHggSelector.hh"


SusyHggSelector::SusyHggSelector(std::vector<std::string> fNames, std::string treeName,std::string outputFile):
  BaseSelector(fNames,treeName,outputFile) 
{
  setDoFill(false);
}

void SusyHggSelector::processConfig(ReadConfig& cfg) {
}

int SusyHggSelector::init() {
  if(optimize) {
    min_pho1_pt=0;
    min_pho2_pt=0;
  }

  if(isMC) {
    jecReader.setCorrections("/home/amott/HggApp/JEC/Summer13_V5_Uncertainties/Summer13_V5_MC_Uncertainty_AK5PFchs.txt");
  }else{
    jecReader.setCorrections("/home/amott/HggApp/JEC/Summer13_V5_Uncertainties/Summer13_V5_DATA_Uncertainty_AK5PFchs.txt");
  }
  assert(jecReader.isValid());
return 0;
}

void SusyHggSelector::clear() {
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

  ele1_pt=0;
  mu1_pt=0;

  nJ=0;
  MR=0;
  Rsq=0;

  nSusyPart=0;
  for(int i=0;i<20;i++) {
    mSusyPart[i]=0;
    idSusyPart[i]=0;
  }

  m22=0;
  m23=0;
  m24=0;
  m25=0;
}


void SusyHggSelector::processEntry(Long64_t iEntry) {
  if(nVtx==0) return;
  
  TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);

  puWeight=pileupWeight;
  runNum  = runNumber;

  //identify the two highest pT photons
  const int NPho=2;
  float pho_sum_pt=0;
  float M_max=-1;

  int selected_photons[NPho];
  for(int i=0;i<NPho;i++) {
    selected_photons[i]=0;
  }
  for(int iPho=0; iPho < nPho_;iPho++) {
    auto photon1 = Photons_->at(iPho);
    //veto photons outside acceptance
    if(fabs(photon1.SC.eta) > max_pho_eta) continue;
    if(fabs(photon1.SC.eta) > 1.4442 && fabs(photon1.SC.eta) < 1.566) continue; //gap photons

    float pt1 = photon1.energy/cosh(photon1.eta);
    if(!optimize) {
      if(pt1 < min_pho2_pt) continue; //not hard enough
      if( !photonID.passID(photon1,StandardPhotonID::kLoose) ) continue; //no photon ID

      std::bitset<5> id_res = photonID.cutResults(photon1,StandardPhotonID::kLoose);
      //if(id_res[3] || id_res[2]) continue; //fails the charged hadron or neutral hadron ID
      
      if(photonMatchedElectron[iPho]) continue; //conversion-safe electron veto
    }
    
    TLorentzVector p4_pho1 = photon1.p4FromVtx(vtx,photon1.finalEnergy);

    for(int jPho=iPho+1;jPho<nPho_;jPho++) {//second photon
      auto photon2 = Photons_->at(jPho);
      //veto photons outside acceptance
      if(fabs(photon2.SC.eta) > max_pho_eta) continue;
      if(fabs(photon2.SC.eta) > 1.4442 && fabs(photon1.SC.eta) < 1.566) continue; //gap photons
      
      float pt2 = photon2.energy/cosh(photon2.eta);
      if(!optimize) {
	if(pt2 < min_pho2_pt) continue; //not hard enough
	if(pt1 <min_pho1_pt && pt2 < min_pho1_pt) continue; //neither is hard enough

	if( !photonID.passID(photon2,StandardPhotonID::kLoose) ) continue; //no photon ID	

	std::bitset<5> id_res = photonID.cutResults(photon2,StandardPhotonID::kLoose);
	//if(id_res[3] || id_res[2]) continue; //fails the charged hadron or neutral hadron ID
	
	if(photonMatchedElectron[jPho]) continue; //conversion-safe electron veto
      }
      TLorentzVector p4_pho2 = photon2.p4FromVtx(vtx,photon2.finalEnergy);
      float M = (p4_pho1+p4_pho2).M();
      if( M > min_mgg && M < max_mgg && (p4_pho1.Pt() + p4_pho2.Pt()) > pho_sum_pt) {
	selected_photons[0]=(pt1>pt2 ? iPho : jPho);
	selected_photons[1]=(pt1>pt2 ? jPho : iPho);
	pho_sum_pt = (p4_pho1.Pt() + p4_pho2.Pt());
      }
    }
  }

  //std::cout << "pho pts: " << pho_pts[0] << "  "<< pho_pts[1] << std::endl;
  //pho_pts is now the sorted pts of the highest pt photons in the event

  //selected_photons[] now contains the (sorted) reference to the two highest pT photons

  TLorentzVector pho1_p4 = Photons_->at(selected_photons[0]).p4FromVtx(vtx,Photons_->at(selected_photons[0]).finalEnergy);
  TLorentzVector pho2_p4 = Photons_->at(selected_photons[1]).p4FromVtx(vtx,Photons_->at(selected_photons[1]).finalEnergy);

  TLorentzVector gg_p4 = pho1_p4+pho2_p4;

  mgg = gg_p4.M();
  //std::cout << "mgg: " << mgg << std::endl;
  ptgg = gg_p4.Pt();
  if(!optimize) {
    if(mgg < min_mgg) return; //mgg too low
    if(ptgg < min_ptgg) return;
  }
  etagg = gg_p4.Eta();
  phigg = gg_p4.Phi();
  pho1_pt = pho1_p4.Pt();
  pho1_eta = pho1_p4.Eta();
  pho1_phi = pho1_p4.Phi();
  pho1_r9 = Photons_->at(selected_photons[0]).SC.r9;
  if(isMC)pho1_seoe = Photons_->at(selected_photons[0]).correctedEnergyError/Photons_->at(selected_photons[0]).correctedEnergy;
  else pho1_seoe = Photons_->at(selected_photons[0]).finalEnergyError/Photons_->at(selected_photons[0]).finalEnergy;

  std::bitset<5> id_res_pho1 = photonID.cutResults(Photons_->at(selected_photons[0]),StandardPhotonID::kLoose);

  pho1_pass_id  = !(id_res_pho1[0] || id_res_pho1[1]);
  pho1_pass_iso = !(id_res_pho1[2] || id_res_pho1[3] || id_res_pho1[4]);


  pho1_sieie   = Photons_->at(selected_photons[0]).SC.sigmaIEtaIEta;
  pho1_HE      = Photons_->at(selected_photons[0]).HoverE;
  pho1_charged = Photons_->at(selected_photons[0]).dr03ChargedHadronPFIso[0];
  pho1_neutral = Photons_->at(selected_photons[0]).dr03NeutralHadronPFIso;
  pho1_photon  = Photons_->at(selected_photons[0]).dr03PhotonPFIso;
  pho1_eleveto = photonMatchedElectron[selected_photons[0]];
  pho1_genMatch = (Photons_->at(selected_photons[0]).genMatch.index!=-1);


  pho2_pt = pho2_p4.Pt();
  pho2_eta = pho2_p4.Eta();
  pho2_phi = pho2_p4.Phi();
  pho2_r9 = Photons_->at(selected_photons[1]).SC.r9;
  if(isMC) pho2_seoe = Photons_->at(selected_photons[1]).correctedEnergyError/Photons_->at(selected_photons[1]).correctedEnergy;
  else pho2_seoe = Photons_->at(selected_photons[1]).finalEnergyError/Photons_->at(selected_photons[1]).finalEnergy;

  std::bitset<5> id_res_pho2 = photonID.cutResults(Photons_->at(selected_photons[1]),StandardPhotonID::kLoose);

  pho2_pass_id  = !(id_res_pho2[0] || id_res_pho2[1]);
  pho2_pass_iso = !(id_res_pho2[2] || id_res_pho2[3] || id_res_pho2[4]);

  pho2_sieie   = Photons_->at(selected_photons[1]).SC.sigmaIEtaIEta;
  pho2_HE      = Photons_->at(selected_photons[1]).HoverE;
  pho2_charged = Photons_->at(selected_photons[1]).dr03ChargedHadronPFIso[0];
  pho2_neutral = Photons_->at(selected_photons[1]).dr03NeutralHadronPFIso;
  pho2_photon  = Photons_->at(selected_photons[1]).dr03PhotonPFIso;
  pho2_eleveto = photonMatchedElectron[selected_photons[1]];
  pho2_genMatch = (Photons_->at(selected_photons[1]).genMatch.index!=-1);

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

  //loop over the jets
  std::vector<TLorentzVector> selectedJets;
  std::vector<TLorentzVector> selectedJets_up;
  std::vector<TLorentzVector> selectedJets_down;

  selectJets(&selectedJets,0);
  selectJets(&selectedJets_up,1);
  selectJets(&selectedJets_down,-1);

  selectedJets.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
  selectedJets_up.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
  selectedJets_down.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
  if(!optimize) {
    if(selectedJets.size()<2) return;
  }

  TVector3 met;
  met.SetPtEtaPhi(pfMet,0.,pfMetPhi);
  
  MET = pfMet;
  METphi = pfMetPhi;
  nJ = selectedJets.size();
  nJ_up = selectedJets_up.size();
  nJ_down = selectedJets_down.size();
  
  if(selectedJets.size()>=2) {
    TLorentzVector h1(0,0,0,0),h2(0,0,0,0);
    try{
      RazorVariables::CombineJets(selectedJets,h1,h2);
     }catch(TooManyJets& e) {
       std::cout << "TOO MANY JETS IN EVENT" << std::endl;
       return;
     }catch(TooFewJets& e) {
       std::cout << "shouldn't happen......" << std::endl;
       throw e;
    }

    if(h1.Pt()==0) return;

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

    Rsq = MTR/MR;
    Rsq=Rsq*Rsq; //square it!
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

      Rsq_up = MTR/MR_up;
      Rsq_up=Rsq_up*Rsq_up; //square it!
    }else{
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

      Rsq_down = MTR/MR_down;
      Rsq_down=Rsq_down*Rsq_down; //square it!
    }else{
      Rsq_down=-1;
      MR_down=-1;
    }
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


void SusyHggSelector::setupOutputTree() {
  outTree = new TTree("SusyHggTree","");

  outTree->Branch("run",&runNum,"run/I");

  outTree->Branch("mgg",&mgg,"mgg/F");
  outTree->Branch("ptgg",&ptgg);
  outTree->Branch("etagg",&etagg);
  outTree->Branch("phigg",&phigg);

  outTree->Branch("pho1_pt",&pho1_pt);
  outTree->Branch("pho1_eta",&pho1_eta);
  outTree->Branch("pho1_phi",&pho1_phi);
  outTree->Branch("pho1_r9",&pho1_r9);
  outTree->Branch("pho1_sigEoE",&pho1_seoe);
  outTree->Branch("pho1_genMatch",&pho1_genMatch,"pho1_genMatch/B");

  outTree->Branch("pho1_sieie",&pho1_sieie);
  outTree->Branch("pho1_HE",&pho1_HE);
  outTree->Branch("pho1_charged",&pho1_charged);
  outTree->Branch("pho1_neutral",&pho1_neutral);
  outTree->Branch("pho1_photon",&pho1_photon);
  outTree->Branch("pho1_eleveto",&pho1_eleveto,"pho1_eleveto/B");
  outTree->Branch("pho1_pass_id",&pho1_pass_id,"pho1_pass_id/B");
  outTree->Branch("pho1_pass_iso",&pho1_pass_iso,"pho1_pass_iso/B");
  
  
  outTree->Branch("pho2_pt",&pho2_pt);
  outTree->Branch("pho2_eta",&pho2_eta);
  outTree->Branch("pho2_phi",&pho2_phi);
  outTree->Branch("pho2_r9",&pho2_r9);
  outTree->Branch("pho2_sigEoE",&pho2_seoe);
  outTree->Branch("pho2_genMatch",&pho2_genMatch,"pho2_genMatch/B");

  outTree->Branch("pho2_sieie",&pho2_sieie);
  outTree->Branch("pho2_HE",&pho2_HE);
  outTree->Branch("pho2_charged",&pho2_charged);
  outTree->Branch("pho2_neutral",&pho2_neutral);
  outTree->Branch("pho2_photon",&pho2_photon);
  outTree->Branch("pho2_eleveto",&pho2_eleveto,"pho2_eleveto/B");
  outTree->Branch("pho2_pass_id",&pho2_pass_id,"pho2_pass_id/B");
  outTree->Branch("pho2_pass_iso",&pho2_pass_iso,"pho2_pass_iso/B");

  outTree->Branch("ele1_pt",&ele1_pt);
  outTree->Branch("ele1_eta",&ele1_eta);
  outTree->Branch("ele1_phi",&ele1_phi);

  outTree->Branch("mu1_pt",&mu1_pt);
  outTree->Branch("mu1_eta",&mu1_eta);
  outTree->Branch("mu1_phi",&mu1_phi);

  outTree->Branch("highest_csv",&highest_csv);
  outTree->Branch("highest_csv_pt",&highest_csv_pt);
  outTree->Branch("highest_csv_eta",&highest_csv_eta);
  outTree->Branch("highest_csv_phi",&highest_csv_phi);

  outTree->Branch("highest_csv_up",&highest_csv_up);
  outTree->Branch("highest_csv_pt_up",&highest_csv_pt_up);
  outTree->Branch("highest_csv_eta_up",&highest_csv_eta_up);
  outTree->Branch("highest_csv_phi_up",&highest_csv_phi_up);

  outTree->Branch("highest_csv_down",&highest_csv_down);
  outTree->Branch("highest_csv_pt_down",&highest_csv_pt_down);
  outTree->Branch("highest_csv_eta_down",&highest_csv_eta_down);
  outTree->Branch("highest_csv_phi_down",&highest_csv_phi_down);

  outTree->Branch("hem1_pt",&hem1_pt);
  outTree->Branch("hem1_eta",&hem1_eta);
  outTree->Branch("hem1_phi",&hem1_phi);
  outTree->Branch("hem1_M",&hem1_M);

  outTree->Branch("hem1_pt_up",&hem1_pt_up);
  outTree->Branch("hem1_eta_up",&hem1_eta_up);
  outTree->Branch("hem1_phi_up",&hem1_phi_up);
  outTree->Branch("hem1_M_up",&hem1_M_up);

  outTree->Branch("hem1_pt_down",&hem1_pt_down);
  outTree->Branch("hem1_eta_down",&hem1_eta_down);
  outTree->Branch("hem1_phi_down",&hem1_phi_down);
  outTree->Branch("hem1_M_down",&hem1_M_down);

  outTree->Branch("hem2_pt",&hem2_pt);
  outTree->Branch("hem2_eta",&hem2_eta);
  outTree->Branch("hem2_phi",&hem2_phi);
  outTree->Branch("hem2_M",&hem2_M);

  outTree->Branch("hem2_pt_up",&hem2_pt_up);
  outTree->Branch("hem2_eta_up",&hem2_eta_up);
  outTree->Branch("hem2_phi_up",&hem2_phi_up);
  outTree->Branch("hem2_M_up",&hem2_M_up);

  outTree->Branch("hem2_pt_down",&hem2_pt_down);
  outTree->Branch("hem2_eta_down",&hem2_eta_down);
  outTree->Branch("hem2_phi_down",&hem2_phi_down);
  outTree->Branch("hem2_M_down",&hem2_M_down);

  outTree->Branch("MET",&MET);
  outTree->Branch("METPhi",&METphi);

  outTree->Branch("Njets",&nJ,"nJets/I");
  outTree->Branch("MR",&MR,"MR/F");
  outTree->Branch("Rsq",&Rsq);

  outTree->Branch("Njets_up",&nJ_up,"nJets_up/I");
  outTree->Branch("MR_up",&MR_up,"MR_up/F");
  outTree->Branch("Rsq_up",&Rsq_up);

  outTree->Branch("Njets_down",&nJ_down,"nJets_down/I");
  outTree->Branch("MR_down",&MR_down,"MR_down/F");
  outTree->Branch("Rsq_down",&Rsq_down);


  outTree->Branch("pileupWeight",&puWeight);


  outTree->Branch("nSusyPart",&nSusyPart,"nSusyPart/I");
  outTree->Branch("mSusyPart",mSusyPart,"mSusyPart[nSusyPart]/F");
  outTree->Branch("idSusyPart",idSusyPart,"idSusyPart[nSusyPart]/I");

  outTree->Branch("m22",&m22);
  outTree->Branch("m23",&m23);
  outTree->Branch("m24",&m24);
  outTree->Branch("m25",&m25);
}


void SusyHggSelector::fillGenTruth() {
  nSusyPart=0;
  for(int iGen=0;iGen<nGenOthers;iGen++) {
    auto gen = GenOthers->at(iGen);
    if(gen.status != 3) continue;
    if(abs(gen.id) < 1000001 || abs(gen.id) > 1000039) continue; //not a susy particle
    idSusyPart[nSusyPart] = gen.id;
    mSusyPart[nSusyPart]  = gen.mass;
    nSusyPart++;
    if(nSusyPart > maxSusyPart) {
      std::cout << "MORE THAN " << maxSusyPart << " GEN SUSY PARTICLES!" << std::endl;
      break;
    }

    switch(abs(gen.id)) {
    case 1000022:
      m22 = gen.mass;
      break;
    case 1000023:
      m23 = gen.mass;
      break;
    case 1000024:
      m24 = gen.mass;
      break;
    case 1000025:
      m25 = gen.mass;
      break;
    }

  }
}

void SusyHggSelector::selectJets(std::vector<TLorentzVector>* selectedJets, int correction) {
  TLorentzVector pho1_p4,pho2_p4;
  pho1_p4.SetPtEtaPhiM( pho1_pt,pho1_eta,pho1_phi,0);
  pho2_p4.SetPtEtaPhiM( pho2_pt,pho2_eta,pho2_phi,0);

  for(int iJet=0;iJet<nJet_;iJet++) {
    auto jet = Jets_->at(iJet);
    TLorentzVector jet_p4 = jet.getP4();
    if(jet_p4.Pt() < 10 || fabs(jet_p4.Eta()) >5) continue; //basic detector cuts

    float corr = 1;
    if(correction==-1) corr-= jecReader.getCorrection(jet_p4.Pt(),jet_p4.Eta(),JECUReader::kDown);
    if(correction==1) corr+= jecReader.getCorrection(jet_p4.Pt(),jet_p4.Eta(),JECUReader::kUp);

    jet_p4.SetPtEtaPhiM( jet_p4.Pt()*corr, jet_p4.Eta(), jet_p4.Phi(), 0. );

    if(jet_p4.Pt() < min_jet_pt) continue;
    if(fabs(jet_p4.Eta()) > max_jet_eta) continue;
    if(jet_p4.DeltaR(pho1_p4) < 0.5 || jet_p4.DeltaR(pho2_p4) < 0.5) continue; //skip jets matched to the photons


    //simple jet quality cuts
    //if(! jetID.passID(jet,VecbosJetID::kLoose) ) continue; //ignore jet ID
    selectedJets->push_back(jet_p4);

    if(correction==0 && jet.combinedSecondaryVertex > highest_csv) { //store the highest csv btag jet
      highest_csv = jet.combinedSecondaryVertex;
      highest_csv_pt = jet_p4.Pt();
      highest_csv_eta = jet_p4.Eta();
      highest_csv_phi = jet_p4.Phi();
    }

    if(correction==1 && jet.combinedSecondaryVertex > highest_csv_up) { //store the highest csv btag jet
      highest_csv_up = jet.combinedSecondaryVertex;
      highest_csv_pt_up = jet_p4.Pt();
      highest_csv_eta_up = jet_p4.Eta();
      highest_csv_phi_up = jet_p4.Phi();
    }

    if(correction==-1 && jet.combinedSecondaryVertex > highest_csv_down) { //store the highest csv btag jet
      highest_csv_down = jet.combinedSecondaryVertex;
      highest_csv_pt_down = jet_p4.Pt();
      highest_csv_eta_down = jet_p4.Eta();
      highest_csv_phi_down = jet_p4.Phi();
    }
  }

}
