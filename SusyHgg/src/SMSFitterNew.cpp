#include "SMSFitterNew.hh"
#include "TH1D.h"
#include <cmath>
#include <iostream>

SMSFitter::SMSFitter(TString inputFileName,TString outputFileName) :Fitter(inputFileName,outputFileName) 
{
  isSMS=false;
  getSMSPoints(&sms_points);
  buildHistograms();
}

void SMSFitter::Run() {
  Long64_t iEntry=-1;

  while(fChain->GetEntry(++iEntry)) {
    int N_pt = nEntriesHist->GetBinContent( nEntriesHist->FindFixBin(m23,m22) );
    assert(N_pt!=0);
    weight = pileupWeight * target_xsec*lumi*HggBR/N_pt;
    processEntry();
  }

  outputFile->cd();
  for(auto hist: SignalRegionHistograms) hist.second->Write();
  for(auto hist: SignalRegionHistogramsFineBin) hist.second->Write();
  outputFile->Close();
}

void SMSFitter::getSMSPoints(std::vector<TString>* pts) {
    for(int m_h=125;m_h<501;m_h+=25) {
      for(int m_l=0;m_l<m_h-124;m_l+=25){
	pts->push_back(Form("%d_%d",m_l,m_h));
      }         
    }
}


void SMSFitter::setNEntriesFile(TString fileName) {
  if(nEntriesFile) nEntriesFile->Close();
  nEntriesFile = new TFile(fileName);

  assert( nEntriesFile->IsOpen() && "Error opening NEntries file " && fileName.Data());

  nEntriesHist = (TH2F*)nEntriesFile->Get("Ntotal");

  assert( nEntriesHist != 0 && "Error getting histogram from NEntries file");
}

SMSFitter::~SMSFitter() {
  nEntriesFile->Close();
}


void SMSFitter::buildHistograms() {
  for(auto sms_pt: sms_points) {
    //std::cout << sms_pt << std::endl;
    for(auto cat: catNames) {
      SignalRegionHistograms[sms_pt+"_"+cat] = new TH2F("data_"+sms_pt+"_"+cat+"_SignalRegion","",nXbins-1,xBins,nYbins-1,yBins);
      SignalRegionHistogramsFineBin[sms_pt+"_"+cat] = new TH2F("data_"+sms_pt+"_"+cat+"_SignalRegion_FineBin","",250,0,2500,100,0,1);
      for(auto dir: systematicDir) {
	for(auto sys: systematicNames) {
	  SignalRegionHistograms[sms_pt+"_"+cat+"_"+sys+"_"+dir] = new TH2F("data_"+sms_pt+"_"+cat+"_SignalRegion_"+sys+"_"+dir,"",nXbins-1,xBins,nYbins-1,yBins);
	}
      }
      mgg_dists[sms_pt+"_"+cat] = new TH1D(sms_pt+"_"+cat+"_mgg_dist","",3000,110,140);
    }
  }
}

void SMSFitter::processEntry() {
  TString sms_pt = getSMSPoint(m22,m23);
    { //nominal 
    TLorentzVector pho1;
    TLorentzVector pho2;
  
    pho1.SetPtEtaPhiM(pho1_pt,pho1_eta,pho1_phi,0);
    pho2.SetPtEtaPhiM(pho2_pt,pho2_eta,pho2_phi,0);

    float se1=pho1_sigEoE;
    float se2=pho2_sigEoE;
    
    float btag = highest_csv;
    float mbbH = mbb_NearH;
    float mbbZ = mbb_NearZ;

    TString cat = getCategory(pho1,pho2,se1,se2,btag,mbbH,mbbZ);
    float sigRegWidth = 2*sigmaEffectives[cat];

    if(mgg > 125+isSMS- sigRegWidth && mgg < 125+isSMS + sigRegWidth) {
      //std::cout << sms_pt << std::endl;
      SignalRegionHistograms[sms_pt+"_"+cat]->Fill(MR,Rsq,weight);
      SignalRegionHistogramsFineBin[sms_pt+"_"+cat]->Fill(MR,Rsq,weight);
    }

    
    mgg_dists[sms_pt+"_"+cat]->Fill(mgg,weight);
  }

  for(auto dir: systematicDir) {
    for(auto sys: systematicNames) {
      TLorentzVector pho1;
      TLorentzVector pho2;
      
      pho1.SetPtEtaPhiM(pho1_pt,pho1_eta,pho1_phi,0);
      pho2.SetPtEtaPhiM(pho2_pt,pho2_eta,pho2_phi,0);
      
      float se1=pho1_sigEoE;
      float se2=pho2_sigEoE;
      
      float btag = highest_csv;
      float mbbH = mbb_NearH;
      float mbbZ = mbb_NearZ;
 
      float thisMR = MR;
      float thisRsq = Rsq;
     
      float mass = mgg;

      if(sys == "phoE") { //photon energy systematic
	float err1 = getSysErrPho(pho1_eta,pho1_r9);
	float err2 = getSysErrPho(pho2_eta,pho2_r9);
	if(dir=="Up") {
	  pho1.SetPtEtaPhiM( pho1.Pt()*(1+err1), pho1.Eta(),pho1.Phi(),0);
	  pho2.SetPtEtaPhiM( pho2.Pt()*(1+err2), pho2.Eta(),pho2.Phi(),0);
	  se1*=1/(1+err1);
	  se2*=1/(1+err2);
	}else{
	  pho1.SetPtEtaPhiM( pho1.Pt()*(1-err1), pho1.Eta(),pho1.Phi(),0);
	  pho2.SetPtEtaPhiM( pho2.Pt()*(1-err2), pho2.Eta(),pho2.Phi(),0);	  
	  se1*=1/(1-err1);
	  se2*=1/(1-err2);
	}
     
	mass = (pho1+pho2).M();
      }

      if(sys == "jec") { //jet energy scale changes Rsq and MR
	if(dir=="Up") {
	  thisMR = MR_up;
	  thisRsq = Rsq_up;
	  btag = highest_csv_up; //on the off chance this changes the jet 
	  mbbH = mbb_NearH_up;
	  mbbZ = mbb_NearZ_up;
	}else{
	  thisMR = MR_down;
	  thisRsq = Rsq_down;
	  btag = highest_csv_down; //on the off chance this changes the jet 
	  mbbH = mbb_NearH_down;
	  mbbZ = mbb_NearZ_down;
	}
      }

      if(sys=="btag") {
	float err =0;
	if(highest_csv_pt > SFb_error.rbegin()->first) err = SFb_error.rbegin()->second;
	else{
	  auto SFBe = SFb_error.begin();
	  for(;SFBe != SFb_error.end(); SFBe++) {
	    if(SFBe->first > highest_csv_pt) {
	      SFBe--;
	      err = SFBe->second;
	      break;
	    }
	  }
	} 
	if(dir=="Up") {
	  btag+=err;
	}else{
	  btag-=err;
	}
      }


      TString cat = getCategory(pho1,pho2,se1,se2,btag,mbbH,mbbZ);
      float sigRegWidth = 2*sigmaEffectives[cat];

      if(mass > 125+isSMS- sigRegWidth && mass < 125+isSMS + sigRegWidth) {
	SignalRegionHistograms[ sms_pt+"_"+cat+"_"+sys+"_"+dir ]->Fill(thisMR,thisRsq,weight);	
      }
    }
  }
}


TString SMSFitter::getSMSPoint(float M22,float M23) {

  int m_lsp = round(M22/25.)*25;
  int m_half = round(M23/25.)*25;
  
  return Form("%d_%d",m_lsp,m_half);
}
