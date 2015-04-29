#include "SMSFitterNew.hh"
#include "TH1D.h"
#include <cmath>
#include <iostream>

SMSFitter::SMSFitter(TString inputFileName,TString outputFileName,bool useHT) :Fitter(inputFileName,outputFileName,useHT) 
{
  isFullSIM = (inputFileName.Contains("FSIM")==0);
  isSMS=false;
  getSMSPoints(&sms_points,isFullSIM);
  buildHistograms();
}

void SMSFitter::Run() {
  Long64_t iEntry=-1;

  while(fChain->GetEntry(++iEntry)) {
    int N_pt = nEntriesHist->GetBinContent( nEntriesHist->FindFixBin(m23,m22) );
    assert(N_pt!=0);
    weight = pileupWeight * target_xsec*lumi*HggBR/N_pt*getBTagSF()*triggerEff;
    if(!passBasicSelection()) continue;
    if(!isFullSIM) scalePhotons();
    processEntry();
  }

  outputFile->cd();
  for(auto hist: SignalRegionHistograms) hist.second->Write();
  for(auto hist: SignalRegionHistogramsFineBin) hist.second->Write();
  outputFile->Close();
}

void SMSFitter::getSMSPoints(std::vector<TString>* pts,bool fullSIM) {
#if 0
    for(int m_h=125;m_h<276;m_h+=25) {
      for(int m_l=0;m_l<1;m_l+=25){
	pts->push_back(Form("%d_%d",m_l,m_h));
      }         
    }
#endif

    if(fullSIM) {
      for(int m_h=125;m_h<201;m_h+=25) {
	for(int m_l=0;m_l<m_h-124;m_l+=25){
	  pts->push_back(Form("%d_%d",m_l,m_h));
	}         
      }
    } else { //FastSIM
      for(int m_h=125;m_h<501;m_h+=25) {
	for(int m_l=0;m_l<m_h-124;m_l+=25){
	  pts->push_back(Form("%d_%d",m_l,m_h));
	}         
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
      mgg_dists[sms_pt+"_"+cat] = new TH1D(sms_pt+"_"+cat+"_mgg_dist","",3000,minMgg,maxMgg);
    }
  }
}

void SMSFitter::processEntry(bool doSyst) {
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

    TString cat = getCategory(pho1,pho2,se1,se2,btag,mbbH,mbbZ,pho1_r9,pho2_r9);
    float sigRegWidth = nSigEffSignalRegion*sigmaEffectives[cat];
    float thisMR = MR;
    float thisRsq=Rsq;
    float thisWeight=weight;
    if(useHT) {
      thisMR = HT;
      thisRsq = MET;
    }

    if(mgg > 125 - sigRegWidth && mgg < 126 + sigRegWidth) {
      //std::cout << sms_pt << std::endl;
      SignalRegionHistograms[sms_pt+"_"+cat]->Fill(thisMR,thisRsq,thisWeight);
      SignalRegionHistogramsFineBin[sms_pt+"_"+cat]->Fill(thisMR,thisRsq,thisWeight);
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
      float thisWeight=weight;
      if(useHT) {
	thisMR = HT;
	thisRsq = MET;
      }
     
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
	  if(useHT) {
	    thisMR = HT_up;
	    thisRsq = MET;
	  }

	  btag = highest_csv_up; //on the off chance this changes the jet 
	  mbbH = mbb_NearH_up;
	  mbbZ = mbb_NearZ_up;
	}else{
	  thisMR = MR_down;
	  thisRsq = Rsq_down;
	  if(useHT) {
	    thisMR = HT_down;
	    thisRsq = MET;
	  }
	  btag = highest_csv_down; //on the off chance this changes the jet 
	  mbbH = mbb_NearH_down;
	  mbbZ = mbb_NearZ_down;
	}
      }

      if(sys=="btag") {
	float btagWeight = getBTagSF();
	if(btagWeight!=1) {
	  float ErrHighest=0;
	  float ErrSecond=0;
	  int iHighest=-1,iSecond=-1;
	  for(int i=0;i<SFbErr_ptMax.size();i++){
	    if(iHighest==-1 && highest_csv_pt < SFbErr_ptMax.at(i)) iHighest=i;
	    if(iSecond==-1 && second_csv_pt < SFbErr_ptMax.at(i)) iSecond=i;	  
	  }
	  //error for the highest CSV 
	  if(highest_csv>0.679) {
	    if(iHighest==-1) {
	      ErrHighest = SFbErr_CSVM.at(SFbErr_CSVM.size()-1)*2;
	    }else{
	      ErrHighest = SFbErr_CSVM.at(iHighest);
	    }	 
	  }else{
	    if(iHighest==-1) {
	      ErrHighest = SFbErr_CSVL.at(SFbErr_CSVL.size()-1)*2;
	    }else{
	      ErrHighest = SFbErr_CSVL.at(iHighest);
	    }	 
	  }

	  //error for the second CSV 
	  if(second_csv>0.679) {
	    if(iSecond==-1) {
	      ErrSecond = SFbErr_CSVM.at(SFbErr_CSVM.size()-1)*2;
	    }else{
	      ErrSecond = SFbErr_CSVM.at(iSecond);
	    }	 
	  }else{
	    if(iSecond==-1) {
	      ErrSecond = SFbErr_CSVL.at(SFbErr_CSVL.size()-1)*2;
	    }else{
	      ErrSecond = SFbErr_CSVL.at(iSecond);
	    }	 
	  }
	
	  if(dir=="Up") {
	    thisWeight*=getBTagSF(ErrHighest,ErrSecond)/btagWeight;
	  }else{
	    thisWeight*=getBTagSF(-1*ErrHighest,-1*ErrSecond)/btagWeight;
	  }
	}
      }

      TString cat = getCategory(pho1,pho2,se1,se2,btag,mbbH,mbbZ,pho1_r9,pho2_r9);
      float sigRegWidth = nSigEffSignalRegion*sigmaEffectives[cat];

      if(mass > 125- sigRegWidth && mass < 126 + sigRegWidth) {
	SignalRegionHistograms[ sms_pt+"_"+cat+"_"+sys+"_"+dir ]->Fill(thisMR,thisRsq,thisWeight);	
      }
    }
  }
}

void SMSFitter::scalePhotons() { //observe an energy shift in the endcap photons -->rescale them
  float p1_scale = getPhotonScale(pho1_eta,pho1_r9);
  pho1_pt*=p1_scale;
  pho1_sigEoE/=p1_scale;

  float p2_scale = getPhotonScale(pho2_eta,pho2_r9);
  pho2_pt*=p2_scale;
  pho2_sigEoE/=p2_scale;
}
float SMSFitter::getPhotonScale(float eta,float r9) {
  if(fabs(eta)< 1.56) return 1;
  if(fabs(eta)< 2.) { //low EE
    if(r9>0.94) return 1.0276;
    else return 1.0216;
  }
  //high EE
  if(r9>0.94) return 1.0287;
  return 1.0351;

}


TString SMSFitter::getSMSPoint(float M22,float M23) {

  int m_lsp = round(M22/25.)*25;
  int m_half = round(M23/25.)*25;
  
  return Form("%d_%d",m_lsp,m_half);
}
