#include "SusyHggSelector.hh"

#include "HggSelector.hh"

SusyHggSelector::SusyHggSelector(std::vector<std::string> fNames, std::string treeName,std::string outputFile):
  BaseSelector(fNames,treeName,outputFile) 
{
  setDoFill(false);
}

void SusyHggSelector::processConfig(ReadConfig& cfg) {
}

void SusyHggSelector::setSmearingCFG(std::string smearCFG) {
  if(smearer!=0) delete smearer;
  smearer = new HggEnergyScale(smearCFG);
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

  if(isDY) {
    min_mgg = 70;
  }
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
  t1Rsq=0;

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


void SusyHggSelector::processEntry(Long64_t iEntry) {
  if(nVtx==0) return;
  
  TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);

  puWeight=pileupWeight;
  runNum  = runNumber;
  lumiSec = lumiBlock;
  evtNum  = evtNumber;
    

    for(int i=0;i<nPho_;i++) {
      auto photon = &(Photons_->at(i));
      if(isMC && smearer) { //smear the photons
	photon->scaledEnergy = photon->correctedEnergy;
	photon->scaledEnergyError = photon->correctedEnergyError;
	std::pair<float,float> dE = smearer->getDEoE(*photon,0);
	photon->dEoE = dE.first;
	photon->dEoEErr = dE.second;

	assert(photon!=0);
	HggSelector::smearPhoton(photon,0);
	//std::cout <<">>>> "<< photon->correctedEnergy << "  " << photon->finalEnergy << "  " << photon->correctedEnergyError << "  " << photon->finalEnergyError  << "  " << photon->dEoE << std::endl;
      }else{
	photon->finalEnergy = photon->correctedEnergy;
	photon->finalEnergyError = photon->correctedEnergyError;
      }
  }

  //identify the two highest pT photons
  const int NPho=2;
  float pho_sum_pt=0;
  float M_max=-1;

  int selected_photons[NPho];
  for(int i=0;i<NPho;i++) {
    selected_photons[i]=0;
  }
  for(int iPho=0; iPho < nPho_;iPho++) {
    auto photon1 =&(Photons_->at(iPho));
    //std::cout << photon1->correctedEnergy << "  " << photon1->finalEnergy << "  " << photon1->correctedEnergyError << "  " << photon1->finalEnergyError  << "  " << photon1->dEoE << std::endl;
    //veto photons outside acceptance
    if(fabs(photon1->SC.eta) > max_pho_eta) continue;
    if(fabs(photon1->SC.eta) > 1.4442 && fabs(photon1->SC.eta) < 1.566) continue; //gap photons

    float pt1 = photon1->energy/cosh(photon1->eta);
    if(!optimize) {
      if(pt1 < min_pho2_pt) continue; //not hard enough
      if( !photonID.passID(*photon1,StandardPhotonID::kLoose) ) continue; //no photon ID

      std::bitset<5> id_res = photonID.cutResults(*photon1,StandardPhotonID::kLoose);
      //if(id_res[3] || id_res[2]) continue; //fails the charged hadron or neutral hadron ID
      
      if(!isDY && photonMatchedElectron[iPho]) continue; //conversion-safe electron veto
    }
    
    TLorentzVector p4_pho1 = photon1->p4FromVtx(vtx,photon1->finalEnergy);

    for(int jPho=iPho+1;jPho<nPho_;jPho++) {//second photon
      auto photon2 = &(Photons_->at(jPho));
      //veto photons outside acceptance
      if(fabs(photon2->SC.eta) > max_pho_eta) continue;
      if(fabs(photon2->SC.eta) > 1.4442 && fabs(photon1->SC.eta) < 1.566) continue; //gap photons
      
      float pt2 = photon2->energy/cosh(photon2->eta);
      if(!optimize) {
	if(pt2 < min_pho2_pt) continue; //not hard enough
	if(pt1 <min_pho1_pt && pt2 < min_pho1_pt) continue; //neither is hard enough

	if( !photonID.passID(*photon2,StandardPhotonID::kLoose) ) continue; //no photon ID	

	std::bitset<5> id_res = photonID.cutResults(*photon2,StandardPhotonID::kLoose);
	//if(id_res[3] || id_res[2]) continue; //fails the charged hadron or neutral hadron ID
	
	if(photonMatchedElectron[jPho]) continue; //conversion-safe electron veto
      }
      TLorentzVector p4_pho2 = photon2->p4FromVtx(vtx,photon2->finalEnergy);
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
  //if(isMC)pho1_seoe = Photons_->at(selected_photons[0]).correctedEnergyError/Photons_->at(selected_photons[0]).correctedEnergy;
  pho1_seoe = Photons_->at(selected_photons[0]).finalEnergyError/Photons_->at(selected_photons[0]).finalEnergy;

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

  pho1_energyGen = (pho1_genMatch ? Photons_->at(selected_photons[0]).genMatch.energy : -1);

  pho2_pt = pho2_p4.Pt();
  pho2_eta = pho2_p4.Eta();
  pho2_phi = pho2_p4.Phi();
  pho2_r9 = Photons_->at(selected_photons[1]).SC.r9;
  //if(isMC) pho2_seoe = Photons_->at(selected_photons[1]).correctedEnergyError/Photons_->at(selected_photons[1]).correctedEnergy;
  pho2_seoe = Photons_->at(selected_photons[1]).finalEnergyError/Photons_->at(selected_photons[1]).finalEnergy;

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

  pho2_energyGen = (pho2_genMatch ? Photons_->at(selected_photons[1]).genMatch.energy : -1);

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

  selectedJets.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
  selectedJets_up.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
  selectedJets_down.push_back(gg_p4); //ensure the 'Higgs' vector is in there explicitly
  if(!optimize) {
    if(selectedJets.size()<2) return;
  }

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
  outTree->Branch("lumi",&lumiSec,"lumi/I");
  outTree->Branch("evt",&evtNum,"evt/I");

  outTree->Branch("mgg",&mgg,"mgg/F");
  outTree->Branch("ptgg",&ptgg);
  outTree->Branch("etagg",&etagg);
  outTree->Branch("phigg",&phigg);
  outTree->Branch("hemgg",&hemgg,"hemgg/I");


  outTree->Branch("pho1_pt",&pho1_pt,"pho1_pt/F");
  outTree->Branch("pho1_eta",&pho1_eta);
  outTree->Branch("pho1_phi",&pho1_phi);
  outTree->Branch("pho1_r9",&pho1_r9);
  outTree->Branch("pho1_sigEoE",&pho1_seoe);
  outTree->Branch("pho1_energyGen",&pho1_energyGen);
  outTree->Branch("pho1_genMatch",&pho1_genMatch,"pho1_genMatch/B");

  outTree->Branch("pho1_sieie",&pho1_sieie,"pho1_sieie/F");
  outTree->Branch("pho1_HE",&pho1_HE);
  outTree->Branch("pho1_charged",&pho1_charged);
  outTree->Branch("pho1_neutral",&pho1_neutral);
  outTree->Branch("pho1_photon",&pho1_photon);
  outTree->Branch("pho1_eleveto",&pho1_eleveto,"pho1_eleveto/B");
  outTree->Branch("pho1_pass_id",&pho1_pass_id,"pho1_pass_id/B");
  outTree->Branch("pho1_pass_iso",&pho1_pass_iso,"pho1_pass_iso/B");
  
  
  outTree->Branch("pho2_pt",&pho2_pt,"pho2_pt/F");
  outTree->Branch("pho2_eta",&pho2_eta);
  outTree->Branch("pho2_phi",&pho2_phi);
  outTree->Branch("pho2_r9",&pho2_r9);
  outTree->Branch("pho2_sigEoE",&pho2_seoe);
  outTree->Branch("pho2_energyGen",&pho2_energyGen);
  outTree->Branch("pho2_genMatch",&pho2_genMatch,"pho2_genMatch/B");

  outTree->Branch("pho2_sieie",&pho2_sieie,"pho2_sieie/F");
  outTree->Branch("pho2_HE",&pho2_HE);
  outTree->Branch("pho2_charged",&pho2_charged);
  outTree->Branch("pho2_neutral",&pho2_neutral);
  outTree->Branch("pho2_photon",&pho2_photon);
  outTree->Branch("pho2_eleveto",&pho2_eleveto,"pho2_eleveto/B");
  outTree->Branch("pho2_pass_id",&pho2_pass_id,"pho2_pass_id/B");
  outTree->Branch("pho2_pass_iso",&pho2_pass_iso,"pho2_pass_iso/B");

  outTree->Branch("ele1_pt",&ele1_pt,"ele1_pt/F");
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

  outTree->Branch("second_csv",&second_csv);
  outTree->Branch("second_csv_pt",&second_csv_pt);
  outTree->Branch("second_csv_eta",&second_csv_eta);
  outTree->Branch("second_csv_phi",&second_csv_phi);

  outTree->Branch("second_csv_up",&second_csv_up);
  outTree->Branch("second_csv_pt_up",&second_csv_pt_up);
  outTree->Branch("second_csv_eta_up",&second_csv_eta_up);
  outTree->Branch("second_csv_phi_up",&second_csv_phi_up);

  outTree->Branch("second_csv_down",&second_csv_down);
  outTree->Branch("second_csv_pt_down",&second_csv_pt_down);
  outTree->Branch("second_csv_eta_down",&second_csv_eta_down);
  outTree->Branch("second_csv_phi_down",&second_csv_phi_down);

  outTree->Branch("mbb",&mbb);
  outTree->Branch("mbb_up",&mbb_up);
  outTree->Branch("mbb_down",&mbb_down);

  outTree->Branch("mbb_NearH",&mbb_NearH);
  outTree->Branch("mbb_NearH_up",&mbb_NearH_up);
  outTree->Branch("mbb_NearH_down",&mbb_NearH_down);

  outTree->Branch("mbb_NearZ",&mbb_NearZ);
  outTree->Branch("mbb_NearZ_up",&mbb_NearZ_up);
  outTree->Branch("mbb_NearZ_down",&mbb_NearZ_down);

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

  outTree->Branch("t1MET",&t1MET);
  outTree->Branch("t1METPhi",&t1METphi);

  outTree->Branch("MR",&MR,"MR/F");
  outTree->Branch("Rsq",&Rsq);
  outTree->Branch("t1Rsq",&t1Rsq);

  outTree->Branch("Njets_up",&nJ_up,"nJets_up/I");
  outTree->Branch("MR_up",&MR_up,"MR_up/F");
  outTree->Branch("Rsq_up",&Rsq_up);
  outTree->Branch("t1Rsq_up",&t1Rsq_up);

  outTree->Branch("Njets_down",&nJ_down,"nJets_down/I");
  outTree->Branch("MR_down",&MR_down,"MR_down/F");
  outTree->Branch("Rsq_down",&Rsq_down);
  outTree->Branch("t1Rsq_down",&t1Rsq_down);


  outTree->Branch("pileupWeight",&puWeight);

  outTree->Branch("HT",&HT,"HT/F");
  outTree->Branch("HT_up",&HT_up,"HT_up/F");
  outTree->Branch("HT_down",&HT_down,"HT_down/F");

  outTree->Branch("MHT",&MHT,"MHT/F");
  outTree->Branch("MHT_up",&MHT_up,"MHT_up/F");
  outTree->Branch("MHT_down",&MHT_down,"MHT_down/F");

  outTree->Branch("MHTphi",&MHTphi,"MHTphi/F");
  outTree->Branch("MHTphi_up",&MHTphi_up,"MHTphi_up/F");
  outTree->Branch("MHTphi_down",&MHTphi_down,"MHTphi_down/F");

  outTree->Branch("Njets",&nJ,"Njets/I");
  outTree->Branch("indexJet",indexJet,"indexJet[Njets]/I");
  outTree->Branch("ptJet",ptJet,"ptJet[Njets]/F");
  outTree->Branch("etaJet",etaJet,"etaJet[Njets]/F");
  outTree->Branch("phiJet",phiJet,"phiJet[Njets]/F");
  outTree->Branch("energyJet",energyJet,"energyJet[Njets]/F");
  outTree->Branch("corrUpJet",corrUpJet,"corrUpJet[Njets]/F");
  outTree->Branch("corrDownJet",corrDownJet,"corrDownJet[Njets]/F");
  outTree->Branch("hemJet",hemJet,"hemJet[Njets]/I");


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
  std::vector<float> btags;


  TLorentzVector pho1_p4,pho2_p4;
  pho1_p4.SetPtEtaPhiM( pho1_pt,pho1_eta,pho1_phi,0);
  pho2_p4.SetPtEtaPhiM( pho2_pt,pho2_eta,pho2_phi,0);


    float *csv_highest_ptr=0,*pt_highest_ptr=0,*eta_highest_ptr=0,*phi_highest_ptr=0;
    float *csv_second_ptr=0,*pt_second_ptr=0,*eta_second_ptr=0,*phi_second_ptr=0;
    float *mbb_ptr;

    float *mbb_NearZ_ptr,*mbb_NearH_ptr;

    switch(correction){
    case 0:
      csv_highest_ptr=&highest_csv;
      pt_highest_ptr =&highest_csv_pt;
      eta_highest_ptr=&highest_csv_eta;
      phi_highest_ptr=&highest_csv_phi;

      csv_second_ptr=&second_csv;
      pt_second_ptr =&second_csv_pt;
      eta_second_ptr=&second_csv_eta;
      phi_second_ptr=&second_csv_phi;

      mbb_ptr = &mbb;
      mbb_NearZ_ptr = &mbb_NearZ;
      mbb_NearH_ptr = &mbb_NearH;
      break;
    case 1:
      csv_highest_ptr=&highest_csv_up;
      pt_highest_ptr =&highest_csv_pt_up;
      eta_highest_ptr=&highest_csv_eta_up;
      phi_highest_ptr=&highest_csv_phi_up;

      csv_second_ptr=&second_csv_up;
      pt_second_ptr =&second_csv_pt_up;
      eta_second_ptr=&second_csv_eta_up;
      phi_second_ptr=&second_csv_phi_up;

      mbb_ptr = &mbb_up;
      mbb_NearZ_ptr = &mbb_NearZ_up;
      mbb_NearH_ptr = &mbb_NearH_up;
      break;
    case -1:
      csv_highest_ptr=&highest_csv_down;
      pt_highest_ptr =&highest_csv_pt_down;
      eta_highest_ptr=&highest_csv_eta_down;
      phi_highest_ptr=&highest_csv_phi_down;

      csv_second_ptr=&second_csv_down;
      pt_second_ptr =&second_csv_pt_down;
      eta_second_ptr=&second_csv_eta_down;
      phi_second_ptr=&second_csv_phi_down;

      mbb_ptr = &mbb_down;
      mbb_NearZ_ptr = &mbb_NearZ_down;
      mbb_NearH_ptr = &mbb_NearH_down;
      break;
    default:
      assert(false);
    }

    

  for(int iJet=0;iJet<nJet_;iJet++) {
    auto jet = Jets_->at(iJet);
    TLorentzVector jet_p4 = jet.getP4();
    if(jet_p4.Pt() < 10 || fabs(jet_p4.Eta()) >5) continue; //basic detector cuts

    float corr = 1;
    try{
      if(correction==-1) corr-= jecReader.getCorrection(jet_p4.Pt(),jet_p4.Eta(),JECUReader::kDown);
      if(correction==1) corr+= jecReader.getCorrection(jet_p4.Pt(),jet_p4.Eta(),JECUReader::kUp);
    }catch(std::runtime_error& e){
      corr=1;
    }
    jet_p4.SetPtEtaPhiM( jet_p4.Pt()*corr, jet_p4.Eta(), jet_p4.Phi(), 0. );

    if(jet_p4.Pt() < min_jet_pt) continue;
    if(fabs(jet_p4.Eta()) > max_jet_eta) continue;
    if(jet_p4.DeltaR(pho1_p4) < 0.5 || jet_p4.DeltaR(pho2_p4) < 0.5) continue; //skip jets matched to the photons


    //simple jet quality cuts
    //if(! jetID.passID(jet,VecbosJetID::kLoose) ) continue; //ignore jet ID
    selectedJets->push_back(jet_p4);
    btags.push_back(jet.combinedSecondaryVertex);

    float thisCSV =jet.combinedSecondaryVertex;
    float thisPt  = jet_p4.Pt();
    float thisEta = jet_p4.Eta();
    float thisPhi = jet_p4.Phi();

    if(thisCSV > *csv_highest_ptr) {
      std::swap(thisCSV,*csv_highest_ptr);
      std::swap(thisPt,*pt_highest_ptr);
      std::swap(thisEta,*eta_highest_ptr);
      std::swap(thisPhi,*phi_highest_ptr);
    }

    if(thisCSV > *csv_second_ptr) {
      std::swap(thisCSV,*csv_second_ptr);
      std::swap(thisPt,*pt_second_ptr);
      std::swap(thisEta,*eta_second_ptr);
      std::swap(thisPhi,*phi_second_ptr);
    }
  }

  *mbb_NearH_ptr=-1;
  *mbb_NearZ_ptr=-1;
  

  for(int ij1=0;ij1<selectedJets->size();ij1++) {
    TLorentzVector j1 = selectedJets->at(ij1);
    if(btags.at(ij1) <0.244) continue;
    for(int ij2=ij1+1;ij2<selectedJets->size();ij2++) {
      TLorentzVector j2 = selectedJets->at(ij2); 
      if(btags.at(ij2) <0.244) continue;
      if(btags.at(ij1) <0.679 && btags.at(ij2)<0.679) continue;
      j1.SetE(sqrt(j1.E()*j1.E()+4.2*4.2));
      j2.SetE(sqrt(j2.E()*j2.E()+4.2*4.2));
      float m = (j1+j2).M();

      if( fabs(*mbb_NearH_ptr-125.) > fabs(m-125.) ) *mbb_NearH_ptr = m;
      if( fabs(*mbb_NearZ_ptr-92.) > fabs(m-92.) ) *mbb_NearZ_ptr = m;
   }
  }

  if(*csv_second_ptr <0) *mbb_ptr=-1;
  else {
    TLorentzVector b1,b2;
    b1.SetPtEtaPhiM(*pt_highest_ptr,*eta_highest_ptr,*phi_highest_ptr,4.2);
    b2.SetPtEtaPhiM(*pt_second_ptr,*eta_second_ptr,*phi_second_ptr,4.2);
  
    *mbb_ptr = (b1+b2).M();
  }
}


/*
int SusyHggSelector::getptbin_for_btag(float pt){
  if(pt<30) return 0;
  else if(pt<40) return 1;
  else if(pt<50) return 2;
  else if(pt<60) return 3;
  else if(pt<70) return 4;
  else if(pt<80) return 5;
  else if(pt<100) return 6;
  else if(pt<120) return 7;
  else if(pt<160) return 8;
  else if(pt<210) return 9;
  else if(pt<260) return 10;
  else if(pt<320) return 11;
  else if(pt<400) return 12;
  else if(pt<500) return 13;
  else if(pt<600) return 14;
  else return 15;
  
}

int SusyHggSelector::get_eta_bin_jet(float eta){
  eta = fabs(eta);
  if(eta<0.9) return 0;
  else if(eta<1.2) return 1;
  else if(eta<2.1) return 2;
  else if(eta<2.4) return 3;
  else return -1;
}

void SusyHggSelector::get_SF_btag(const VecbosJet& jet, float &SF, float &SFerr){
  float SFb_error[100] = {
    0.0415694,
    0.023429,
    0.0261074,
    0.0239251,
    0.0232416,
    0.0197251,
    0.0217319,
    0.0198108,
    0.0193,
    0.0276144,
    0.0205839,
    0.026915,
    0.0312739,
    0.0415054,
    0.0740561,
  0.0598311 
  };

  float x = jet.pt; ///the pt of the jet 
  float eta = fabs(jet.eta); ///abs(eta) 

  if(eta>2.4){
    std::cout<<"warning SF_btag_eta>2.4 ?? " << eta << std::endl; 
    exit(1);
  }
  
  if(x<20) x=20; 
  if(x>800) x= 800;

  //figure out the MC true flavor

  auto matchLambda = [](const VecbosJet& j, std::vector<VecbosGen>* v,float & minDR) {
    index = -1;
    for(int i=0;i<v->size();i++) {
      auto gen = v->at(i);
      if( fabs(j.pt-gen.pt)/gen.pt > 0.5 && fabs(j.pt-gen.pt)/gen.pt < 2. ){
	TLorentzVector gp4; gp4.SetPtEtaPhiM(gen.pt,gen.eta,gen.phi,gen.mass);
	float dr = gp4.DeltaR(j.getP4());
	if(dr < minDR) {
	  minDR = dr;
	  index=i;
	}
      }
    }
    return index;
  }
  
  bool isLight=false;

  float smallestDR = 10;
  int index = matchLambda(j,GenOthers,smallestDR);
  float dr;
  matchLambda(j,GenElectrons,dr);
  matchLambda(j,GenMuons,dr);
  matchLambda(j,GenPhotons,dr);
  if(dr < smallestDR) { isLight=true; } //don't care what flavor of lepton/photon it is, it can't be a b jet
  if(index==-1) isLight=true;

  if(!isLight) {
    int flavor = GenOthers->At(index)->id;
    if( abs(flavor) != 5 && abs(flavor) !=4) isLight=true;
  }

  

  if(!isLight){ //for b or c. flavJet[indj] refers to the MC-true flavor of the jet
    SF  = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x)); 
    int ptbin = getptbin_for_btag( jet.pt ); ///--> ptJet[indj] refers to the pt of this jet
    SFerr = SFb_error[ptbin];
    if(x>800 || x<20 ) SFerr *= 2;
    if(abs(flavJet[indj]) == 4 ) SFerr *= 2; 
  }else{ ///SFlight
    float max;
    float min;
    if(eta<=0.8){
      SF = ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
      max = ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
      min = ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
    }
    else if(eta<=1.6){
      SF = ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
      max = ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
      min = ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
    }  else if(eta>1.6 && eta<=2.4){
      SF = ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
      max = ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
      min = ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
    }
    SFerr = fabs(max-SF)>fabs(min-SF)? fabs(max-SF):fabs(min-SF);
    
  }
  
}
void SusyHggSelector::computeBTagSF() {
  float mcTag = 1.;
  float mcNoTag = 1.;
  float dataTag = 1.;
  float dataNoTag = 1.;
  float errTag = 0.;
  float errNoTag = 0.;

  float err1 = 0; 
  float err2 = 0; 
  flaot err3 = 0; 
  float err4 = 0; 
  
  for(int nj=0; nj<int(indseljet.size());nj++){ //Here we loop over all selected jets ( for our case, pt>30, PF loose ID, etc ) 
    int indj = indseljet[nj];
    float csv = btagJet[indj][0]; ////here we get the CSV btag value 
    int partonFlavor = abs(flavJet[indj]);
    float eta = fabs(etaJet[indj]);
    if(eta>2.4) continue;
    if(partonFlavor==0) continue; //for jets with flavor 0, we ignore. 

    int etabin = get_eta_bin_jet(eta);
    float eff; 
    if( partonFlavor==5 ) {
      ///here one need to provide the pt/eta dependent efficiency for b-tag for "b jet"
    }else if( partonFlavor==4){
      ///here one need to provide the pt/eta dependent efficiency for b-tag for "c jet"
    }else{
      ///here one need to provide the pt/eta dependent efficiency for b-tag for "light jet"
    }
    bool istag = csv > 0.679 && fabs(etaJet[indj])<2.4 ;
    float SF = 0.;
    float SFerr = 0.;
    get_SF_btag(indj,SF,SFerr);

    if(istag){
      mcTag *= eff; 
      dataTag *= eff*SF; 

      if(partonFlavor==5 || partonFlavor ==4)  err1 += SFerr/SF; ///correlated for b/c
      else err3 += SFerr/SF; //correlated for light
      
      
    }else{
      mcNoTag *= (1- eff); 
      dataNoTag *= (1- eff*SF); 
      
      if(partonFlavor==5 || partonFlavor ==4 ) err2 += (-eff*SFerr)/(1-eff*SF); /// /correlated for b/c
      else err4 +=  (-eff*SFerr)/(1-eff*SF);  ////correlated for light
      
    }
    
  }
  
  wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag ); 
  wtbtagErr = sqrt( pow(err1+err2,2) + pow( err3 + err4,2)) * wtbtag;  ///un-correlated for b/c and light
}
*/
