#include "FitterNew.hpp"

#include "TObjArray.h"
#include "TMath.h"

#define NUM_CPU 1

const std::vector<TString> Fitter::catNames = { "HighPt","BTag","HighRes","LowRes"};
const std::vector<TString> Fitter::systematicNames = {"jec","phoE","btag","sigE"}; 


Fitter::Fitter(TString inputFileName,TString outputFileName):SusyHggTree(inputFileName) {
  outputFile = new TFile(outputFileName,"RECREATE");

  buildHistograms();  

  setSigEff("HighPt",3);
  setSigEff("BTag",3);
  setSigEff("HighRes",3);
  setSigEff("LowRes",3);
}

Fitter::~Fitter() {
  outputFile->Close();
}

void Fitter::buildHistograms() {
  for(auto cat: catNames) {
    SignalRegionHistograms[cat] = new TH2F("data_"+cat+"_SignalRegion","",nXbins-1,xBins,nYbins-1,yBins);
    SignalRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SignalRegion_FineBin","",250,0,2500,100,0,1);
    for(auto dir: systematicDir) {
      for(auto sys: systematicNames) {
	SignalRegionHistograms[cat+"_"+sys+"_"+dir] = new TH2F("data_"+cat+"_SignalRegion_"+sys+"_"+dir,"",nXbins-1,xBins,nYbins-1,yBins);
      }
    }
    mgg_dists[cat] = new TH1D(cat+"_mgg_dist","",3000,110,140);
  }
}

bool Fitter::passBasicSelection() {
  if(!pho1_pass_iso || !pho2_pass_iso) return false;
  return true;
}

TString Fitter::getCategory(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag) {
  float ptgg = (pho1+pho2).Pt();

  if(ptgg > 110) return catNames[0];

  if(btag > 0.244) return catNames[1];
  
  if(  se1 > ( fabs(pho1.Eta())>1.48 ? 0.024 : 0.015 ) ) return catNames[3]; //pho1 fails se/e cut
  if(  se2 > ( fabs(pho2.Eta())>1.48 ? 0.024 : 0.015 ) ) return catNames[3]; //pho1 fails se/e cut
  
  return catNames[2]; //high res category
}

void Fitter::processEntry() {
  { //nominal 
    TLorentzVector pho1;
    TLorentzVector pho2;
  
    pho1.SetPtEtaPhiM(pho1_pt,pho1_eta,pho1_phi,0);
    pho2.SetPtEtaPhiM(pho2_pt,pho2_eta,pho2_phi,0);

    float se1=pho1_sigEoE;
    float se2=pho2_sigEoE;
    
    float btag = highest_csv;
    TString cat = getCategory(pho1,pho2,se1,se2,btag);
    float sigRegWidth = 2*sigmaEffectives[cat];

    if(mgg > 125+isSMS- sigRegWidth && mgg < 125+isSMS + sigRegWidth) {
      SignalRegionHistograms[cat]->Fill(MR,Rsq,weight);
      SignalRegionHistogramsFineBin[cat]->Fill(MR,Rsq,weight);
    }

    
    mgg_dists[cat]->Fill(mgg,weight);
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
 
      float thisMR = MR;
      float thisRsq = Rsq;
     
      float mass = mgg;

      if(sys == "phoE") { //photon energy systematic
	float err1 = getSysErrPho(pho1_eta,pho1_r9);
	float err2 = getSysErrPho(pho2_eta,pho2_r9);
	assert(err1 > 0. && err1 < 0.1);
	assert(err2 > 0. && err2 < 0.1);
	if(dir=="Up") {
	  pho1.SetPtEtaPhiM( pho1.Pt()*(1+err1), pho1.Eta(),pho1.Phi(),0);
	  pho2.SetPtEtaPhiM( pho2.Pt()*(1+err2), pho2.Eta(),pho2.Phi(),0);
	}else{
	  pho1.SetPtEtaPhiM( pho1.Pt()*(1-err1), pho1.Eta(),pho1.Phi(),0);
	  pho2.SetPtEtaPhiM( pho2.Pt()*(1-err2), pho2.Eta(),pho2.Phi(),0);	  
	}
     
	mass = (pho1+pho2).M();
      }

      if(sys == "jec") { //jet energy scale changes Rsq and MR
	if(dir=="Up") {
	  thisMR = MR_up;
	  thisRsq = Rsq_up;
	  btag = highest_csv_up; //on the off chance this changes the jet 
	}else{
	  thisMR = MR_down;
	  thisRsq = Rsq_down;
	  btag = highest_csv_down; //on the off chance this changes the jet 
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


      TString cat = getCategory(pho1,pho2,se1,se2,btag);
      float sigRegWidth = 2*sigmaEffectives[cat];

      if(mass > 125+isSMS- sigRegWidth && mass < 125+isSMS + sigRegWidth) {
	SignalRegionHistograms[ cat+"_"+sys+"_"+dir ]->Fill(thisMR,thisRsq,weight);	
      }
    }
  }
}

float Fitter::getSysErrPho(float eta,float r9) {
     TString region = "EBLow";
     TString r9s = (r9 > 0.94 ? "highR9" : "lowR9");

     if(fabs(eta) > 1.0)  region = "EBHigh";
     if(fabs(eta) > 1.48) region = "EELow";
     if(fabs(eta) > 2.0)  region = "EEHigh";

     return smearSys[ Form("%s_%s",region.Data(),r9s.Data()) ];
}


void Fitter::Run() {
  Long64_t iEntry=-1;
  while(fChain->GetEntry(++iEntry)) {
    weight = pileupWeight * target_xsec*lumi*HggBR/N_total;
    processEntry();
  }

  outputFile->cd();
  for(auto hist: SignalRegionHistograms) hist.second->Write();
  for(auto hist: SignalRegionHistogramsFineBin) hist.second->Write();
  for(auto hist: mgg_dists) hist.second->Write();
  outputFile->Close();
}

