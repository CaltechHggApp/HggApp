#include "FitterNew.hpp"

#include "TObjArray.h"
#include "TMath.h"

#include <iostream>

#define NUM_CPU 1

//#define USER9 //use r9 base categories

const std::vector<TString> Fitter::catNames = { "HighPt","Hbb","Zbb","HighRes","LowRes"};
const std::vector<TString> Fitter::systematicNames = {"jec","phoE","btag","sigE","statistics"}; 
//const std::vector<TString> Fitter::catNames = { "Hbb"};
//const std::vector<TString> Fitter::systematicNames = {}; 


Fitter::Fitter(TString inputFileName,TString outputFileName,bool useHT):SusyHggTree(inputFileName) {
  for(int i=0;i<nZbins;i++) { zBins[i] = 100+0.5*i; }

  outputFile = new TFile(outputFileName,"RECREATE");
  setUseHT(useHT);
  buildHistograms();  


  for(auto cat: catNames) {
    setSigEff(cat,3);
    nSignal[cat] = 0;
    nTotal[cat] = 0;
  }
}

Fitter::~Fitter() {
  outputFile->Close();
}

void Fitter::buildHistograms() {
  for(auto cat: catNames) {
    if(useHT) {
      SignalRegionHistograms[cat] = new TH2F("data_"+cat+"_SignalRegion","",47,150,2500,12,0,300);
      SignalRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SignalRegion_FineBin","",500,0,2500,100,0,300);
      for(auto dir: systematicDir) {
	for(auto sys: systematicNames) {
	  SignalRegionHistograms[cat+"_"+sys+"_"+dir] = new TH2F("data_"+cat+"_SignalRegion_"+sys+"_"+dir,"",47,150,2500,12,0,300);
	}
      }
      SignalRegions3D[cat] = new TH3F("data_"+cat+"_SignalRegion_3D","",47,150,2500,12,0,300,200,100,200);
    } else {
      SignalRegionHistograms[cat] = new TH2F("data_"+cat+"_SignalRegion","",nXbins-1,xBins,nYbins-1,yBins);
      SignalRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SignalRegion_FineBin","",250,0,2500,100,0,1);
      for(auto dir: systematicDir) {
	for(auto sys: systematicNames) {
	  SignalRegionHistograms[cat+"_"+sys+"_"+dir] = new TH2F("data_"+cat+"_SignalRegion_"+sys+"_"+dir,"",nXbins-1,xBins,nYbins-1,yBins);
	}
      }
      SignalRegions3D[cat] = new TH3F("data_"+cat+"_SignalRegion_3D","",nXbins-1,xBins,nYbins-1,yBins,nZbins-1,zBins);

    }
    mgg_dists[cat] = new TH1D(cat+"_mgg_dist","",3000,minMgg,maxMgg);
  }
}

bool Fitter::passBasicSelection() {
  if(!pho1_pass_iso || !pho2_pass_iso) return false;
  if(ptgg<10) return false;
  if( !alternateAnalysis && ptgg<20 ) return false;

  if(mgg<100 || mgg>180) return false;

  //if(ele1_pt>15 || mu1_pt>15) return false;

  switch(basicSelection) {
  case kAN239:
    //AN13/239-like selection
    if(pho1_pt<40 || pho2_pt < 25) return false;
    if( fabs(pho1_eta)>1.46 || fabs(pho2_eta)>1.46) return false;
    break;

  case kHighPt:
    if(pho1_pt<40 || pho2_pt < 25) return false;
    break;

  case kLoose:
    if(pho1_pt<32 || pho2_pt < 24) return false;
    break;

  default:
    std::cout << "invalid basic selection (or selection not implemented)" << std::endl;
    assert(false);

  }
  return true;
}

TString Fitter::getCategory(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag,float mbbH,float mbbZ,float r9_1,float r9_2) {
  if(alternateAnalysis) return getCategoryAlt(pho1,pho2,se1,se2,btag,mbbH,mbbZ,r9_1,r9_2);
  return getCategoryOrig(pho1,pho2,se1,se2,btag,mbbH,mbbZ,r9_1,r9_2);
}

TString Fitter::getCategoryOrig(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag,float mbbH,float mbbZ,float r9_1,float r9_2) {
  float ptgg = (pho1+pho2).Pt();

  //return catNames[0];

  if(ptgg > 110) return catNames[0];

 
  
  if( fabs(mbbH-125.) <15. ) return catNames[1];
  if( fabs(mbbZ-91.) <15. ) return catNames[2];
  
  if(  se1 > ( fabs(pho1.Eta())>1.48 ? 0.024 : 0.015 ) ) return catNames[4]; //pho1 fails se/e cut
  if(  se2 > ( fabs(pho2.Eta())>1.48 ? 0.024 : 0.015 ) ) return catNames[4]; //pho1 fails se/e cut
  
  return catNames[3]; //high res category

}

TString Fitter::getCategoryAlt(const TLorentzVector& pho1, const TLorentzVector&pho2,float se1, float se2,float btag,float mbbH,float mbbZ,float r9_1,float r9_2) {
  float ptgg = (pho1+pho2).Pt();
  if(ptgg > 110) return catNames[0];
  if( mbbH>108 && mbbH<150 ) return catNames[1];
  if( mbbZ>66  && mbbZ<108 ) return catNames[2];
  
  if(r9_1<0.94 || r9_2<0.94) return catNames[4];

  return catNames[3]; //high res category
}

void Fitter::processEntry(bool doSyst) {
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
    float thisRsq = Rsq;
    float thisWeight=weight;
    assert(!isnan(thisWeight));
    if(useHT) {
      thisMR = HT;
      thisRsq = MET;
    }

    if(mgg > 125- sigRegWidth && mgg < 126 + sigRegWidth) {
      SignalRegionHistograms[cat]->Fill(thisMR,thisRsq,thisWeight);
      SignalRegionHistogramsFineBin[cat]->Fill(thisMR,thisRsq,thisWeight);
      SignalRegions3D[cat]->Fill(thisMR,thisRsq,mgg,thisWeight);
    }
    nSignal[cat]++;
    nTotal[cat]++;
    
    mgg_dists[cat]->Fill(mgg,thisWeight);
  }

  if(!doSyst)  return;

  //do systematics
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
	assert(err1 > 0. && err1 < 0.1);
	assert(err2 > 0. && err2 < 0.1);
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
	assert(btagWeight>0);
	assert(!isnan(btagWeight));
	assert(SFbErr_ptMax.size() == SFbErr_CSVL.size());
	assert(SFbErr_ptMax.size() == SFbErr_CSVM.size());
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
	  if(isnan(ErrHighest) || isnan(ErrSecond)) {
	    std::cerr << "a btag weight is NAN: " << ErrHighest << " " << ErrSecond << std::endl;
	    assert(false);
	  }
	  if(dir=="Up") {
	    thisWeight*=getBTagSF(ErrHighest,ErrSecond)/btagWeight;
	  }else{
	    thisWeight*=getBTagSF(-1*ErrHighest,-1*ErrSecond)/btagWeight;
	  }
	  assert(!isnan(thisWeight));
	}
      }


      TString cat = getCategory(pho1,pho2,se1,se2,btag,mbbH,mbbZ,pho1_r9,pho2_r9);
      float sigRegWidth = nSigEffSignalRegion*sigmaEffectives[cat];

      if(mass > 125- sigRegWidth && mass < 126 + sigRegWidth) {
	SignalRegionHistograms[ cat+"_"+sys+"_"+dir ]->Fill(thisMR,thisRsq,thisWeight);	
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
  TFile* metPhiRatio_f=0;
  TH1D* metPhiRatio=0;
  if(MetPhiSF_file!="") {
    metPhiRatio_f=new TFile(MetPhiSF_file);
    metPhiRatio=(TH1D*)metPhiRatio_f->Get("weight");
  }
  
  while(fChain->GetEntry(++iEntry)) {
    weight = pileupWeight * hggSigStrength*target_xsec*lumi*HggBR*triggerEff/N_total;
    weight*=getBTagSF();
    if(metPhiRatio) {
      weight*= metPhiRatio->GetBinContent(metPhiRatio->FindFixBin(METPhi));
    }
    if(!passBasicSelection()) continue;
    processEntry();
  }

  for(auto cat: catNames) {
    float scale = nSignal[cat]/SignalRegionHistograms[cat]->Integral();
    for(int xBin=0; xBin<SignalRegionHistograms[cat]->GetNbinsX()+1; xBin++) {
      for(int yBin=0; yBin<SignalRegionHistograms[cat]->GetNbinsY()+1; yBin++) {
	float content = SignalRegionHistograms[cat]->GetBinContent(xBin,yBin);
	float error = TMath::Sqrt(content/scale);
	SignalRegionHistograms[ cat+"_statistics_Up" ]->SetBinContent(xBin,yBin,content+error);
	SignalRegionHistograms[ cat+"_statistics_Down" ]->SetBinContent(xBin,yBin,content-error);      
      }
    }
    
  }
  outputFile->cd();
  for(auto hist: SignalRegionHistograms) hist.second->Write();
  for(auto hist: SignalRegionHistogramsFineBin) hist.second->Write();
  for(auto hist: mgg_dists) hist.second->Write();
  for(auto hist: SignalRegions3D) hist.second->Write();
  outputFile->Close();
}

 float Fitter::getBTagSF(float errHigh, float errLow) {
  if(highest_csv < 0.244) return 1;
  if(second_csv<0.244) return 1;
  float sf=1;

  if(highest_csv>0.679) sf*=getSFb(highest_csv_pt,false)+errHigh;
  else sf*=getSFb(highest_csv_pt,true)+errHigh;

  if(second_csv>0.679) sf*=getSFb(second_csv_pt,false)+errLow;
  else sf*=getSFb(second_csv_pt,true)+errLow;
  return sf;
}

float Fitter::getSFb(float pt, bool CSVL) {
  if(pt>800) pt=800;
  if(CSVL) return 0.997942*((1.+(0.00923753*pt))/(1.+(0.0096119*pt)));
  return (0.938887+(0.00017124*pt))+(-2.76366E-07*(pt*pt));
}
