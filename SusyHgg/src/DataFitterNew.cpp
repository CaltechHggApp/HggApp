#include "DataFitterNew.hh"

#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"

#include "assert.h"

RealVar DataFitter::doFitGetScale(TTree* data,float width,RooWorkspace* ws,bool dExp) {

  assert(data->GetEntries() != 0);

  RooRealVar mgg("mgg","",minMgg,maxMgg);

  RooRealVar a1("a1","",-0.3,-2.,0.);
  RooRealVar a2("a2","",-0.1,-2.,0.);
  RooRealVar f("f","",0.5,0.00,1.0);

  RooRealVar NBkg("NBkg","",1000,0,1e9);

  RooExponential e1("e1","",mgg,a1);
  RooExponential e2("e2","",mgg,a2);

  RooAddPdf add("add","",e1,e2,f);

  RooExtendPdf *pdf=0;
  if(dExp) pdf = new RooExtendPdf("pdf","",add,NBkg);
  else     pdf = new RooExtendPdf("pdf","",e1,NBkg);

  RooDataSet rdata("data","",data,mgg);
    

  pdf->fitTo(rdata,RooFit::Strategy(0),RooFit::Extended(kTRUE));
  RooFitResult* res = pdf->fitTo(rdata,RooFit::Strategy(2),RooFit::Extended(kTRUE),RooFit::Save(kTRUE));
  
  mgg.setRange("sig",125-width,126+width);

  RooAbsReal *sig_int = pdf->createIntegral(mgg,RooFit::NormSet(mgg),RooFit::Range("sig"));
  sig_int->SetName("signal_integral");


  float N_sideband = rdata.sumEntries(Form("(mgg>%0.2f && mgg <120) || (mgg>131 && mgg<%0.2f)",minMgg,maxMgg));  

  RealVar scale;
  scale.val = NBkg.getVal()*sig_int->getVal()/N_sideband;
  scale.error = scale.val*sqrt( 1/N_sideband + pow(sig_int->getPropagatedError(*res)/sig_int->getVal(),2) + 1/NBkg.getVal() );

  RooRealVar scaleFactor("scaleFactor","",scale.val);
  scaleFactor.setError(scale.error);

  if(ws) {
    ws->import(*pdf);
    ws->import(rdata);
    ws->import(*res);
    ws->import(*sig_int);
    ws->import(scaleFactor);
  }else{
    delete res;
    delete sig_int;
  }
  return scale;
}

void DataFitter::buildSidebandHistograms() {
  for(auto cat: catNames) {
    if(useHT) {
      SidebandRegionHistograms[cat] = new TH2F("data_"+cat+"_SidebandRegion","",50,0,2500,12,0,300);
      SidebandRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SidebandRegion_FineBin","",500,0,2500,100,0,300);
      SidebandRegionHistograms[cat+"_statistics_Up"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsUp","",50,0,2500,12,0,300);
      SidebandRegionHistograms[cat+"_statistics_Down"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsDown","",50,0,2500,12,0,300);
      SidebandRegionHistograms[cat+"_fit_Up"] = new TH2F("data_"+cat+"_SidebandRegion_fitUp","",50,0,2500,12,0,300);
      SidebandRegionHistograms[cat+"_fit_Down"] = new TH2F("data_"+cat+"_SidebandRegion_fitDown","",50,0,2500,12,0,300);
      SidebandRegionHistograms[cat+"_bkgShape_Up"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeUp","",50,0,2500,12,0,300);
      SidebandRegionHistograms[cat+"_bkgShape_Down"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeDown","",50,0,2500,12,0,300);
    }else{
      SidebandRegionHistograms[cat] = new TH2F("data_"+cat+"_SidebandRegion","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SidebandRegion_FineBin","",250,0,2500,100,0,1);
      SidebandRegionHistograms[cat+"_statistics_Up"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_statistics_Down"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsDown","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_fit_Up"] = new TH2F("data_"+cat+"_SidebandRegion_fitUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_fit_Down"] = new TH2F("data_"+cat+"_SidebandRegion_fitDown","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_bkgShape_Up"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_bkgShape_Down"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeDown","",nXbins-1,xBins,nYbins-1,yBins);
    }

    TTree massTree("massTree","");
    float m;
    massTree.Branch("mgg",&m);
    Long64_t iEntry=-1;
    while( fChain->GetEntry(++iEntry) ) {
      if(!passBasicSelection()) continue;
      if(hasTrigger && !triggerBit) continue;
      TLorentzVector pho1;
      TLorentzVector pho2;
  
      pho1.SetPtEtaPhiM(pho1_pt,pho1_eta,pho1_phi,0);
      pho2.SetPtEtaPhiM(pho2_pt,pho2_eta,pho2_phi,0);
      
      float se1=pho1_sigEoE;
      float se2=pho2_sigEoE;
      
      float btag = highest_csv;

      if( getCategory(pho1,pho2,se1,se2,btag,mbb_NearH,mbb_NearZ,pho1_r9,pho2_r9) == cat ) {
	m=mgg;
	massTree.Fill();
      }
    }
    
    std::cout << cat << std::endl;
    RooWorkspace *ws = new RooWorkspace(cat+"_mgg_workspace","");
    scales[cat] = doFitGetScale(&massTree,nSigEffSignalRegion*sigmaEffectives[cat],ws,(cat=="HighPt" ? false:true) );

    mggFitWorkspaces.push_back(ws);
  }
}

void DataFitter::Run() {
  for( auto cat: catNames) {
    nSideband[cat]=0;
  }

  buildSidebandHistograms();
  Long64_t iEntry=-1;
  Long64_t nPassSelection=0;
  Long64_t nPassTrigger=0;
  Long64_t nPassNoise=0;
  while(fChain->GetEntry(++iEntry)) {
    weight=1; //for data
    if(!passBasicSelection()) continue;
    nPassSelection++;
    if(hasTrigger && !triggerBit) continue;
    nPassTrigger++;
    if(hasNoise && !noiseBit) continue;
    nPassNoise++;
    processEntry();
    processEntrySidebands();
  }
  std::cout << "nProcessed:  " << iEntry << std::endl;
  for(auto cat: catNames) {
    std::cout << cat << "  " << SidebandRegionHistograms[cat]->Integral() << std::endl;

    SidebandRegionHistograms[cat+"_bkgShape_Up"]->Scale( SidebandRegionHistograms[cat]->Integral()/SidebandRegionHistograms[cat+"_bkgShape_Up"]->Integral() );
    SidebandRegionHistograms[cat+"_bkgShape_Down"]->Scale( SidebandRegionHistograms[cat]->Integral()/SidebandRegionHistograms[cat+"_bkgShape_Down"]->Integral() );
  }
  
  for(auto cat: catNames) {
    assert(nSideband[cat]!=0);
    float scale = nSideband[cat]/SidebandRegionHistograms[cat]->Integral(); //scale up to full statistics
    for(int xBin=0; xBin<SidebandRegionHistograms[cat]->GetNbinsX()+1; xBin++) {
      for(int yBin=0; yBin<SidebandRegionHistograms[cat]->GetNbinsY()+1; yBin++) {
	float content = SidebandRegionHistograms[cat]->GetBinContent(xBin,yBin);
	float error = TMath::Sqrt(content/scale);
	SidebandRegionHistograms[ cat+"_statistics_Up" ]->SetBinContent(xBin,yBin,content+error);
	SidebandRegionHistograms[ cat+"_statistics_Down" ]->SetBinContent(xBin,yBin,content-error);      
      }
    }
  }


  outputFile->cd();
  for(auto hist: SignalRegionHistograms) {
    hist.second->Write();
  }
  for(auto hist: SignalRegionHistogramsFineBin) hist.second->Write();

  for(auto hist: SidebandRegionHistograms) hist.second->Write();
  for(auto hist: SidebandRegionHistogramsFineBin) hist.second->Write();

  for(auto ws: mggFitWorkspaces) ws->Write();

  outputFile->Close();

  for(auto scale: scales) std::cout << scale.first << "   " << scale.second.val << " +- " << scale.second.error << std::endl;

  std::cout << nPassSelection << " / " << iEntry << "  events passing baseline selection" << std::endl;
  std::cout << nPassTrigger << " / " << nPassSelection << "  events passing trigger" << std::endl;
  std::cout << nPassNoise << " / " << nPassTrigger << "  events passing noise filters" << std::endl;
}

void DataFitter::processEntrySidebands() {
    TLorentzVector pho1;
    TLorentzVector pho2;
  
    pho1.SetPtEtaPhiM(pho1_pt,pho1_eta,pho1_phi,0);
    pho2.SetPtEtaPhiM(pho2_pt,pho2_eta,pho2_phi,0);

    float se1=pho1_sigEoE;
    float se2=pho2_sigEoE;
    
    float btag = highest_csv;
    TString cat = getCategory(pho1,pho2,se1,se2,btag,mbb_NearH,mbb_NearZ,pho1_r9,pho2_r9);
    float sigRegWidth = nSigEffSignalRegion*sigmaEffectives[cat];

    float thisMR = MR;
    float thisRsq = Rsq;
    if(useHT) {
      thisMR = HT;
      thisRsq = MET;
    }

    if((mgg>minMgg && mgg<120) || (mgg>131 && mgg<maxMgg)){
      SidebandRegionHistograms[cat]->Fill(thisMR,thisRsq,scales[cat].val*weight);
      SidebandRegionHistograms[cat+"_fit_Up"]->Fill(thisMR,thisRsq,(scales[cat].val+scales[cat].error)*weight);
      SidebandRegionHistograms[cat+"_fit_Down"]->Fill(thisMR,thisRsq,(scales[cat].val-scales[cat].error)*weight);
      SidebandRegionHistogramsFineBin[cat]->Fill(thisMR,thisRsq,scales[cat].val*weight);
      nSideband[cat]++;
    }
    if(mgg > 131 && mgg < maxMgg) {
      SidebandRegionHistograms[cat+"_bkgShape_Up"]->Fill(thisMR,thisRsq,scales[cat].val*weight);
    }
    if(mgg > minMgg && mgg < 120) {
      SidebandRegionHistograms[cat+"_bkgShape_Down"]->Fill(thisMR,thisRsq,scales[cat].val*weight);
    }
}

float DataFitter::getSysErrPho(float eta, float r9) {
  TString region = "EBLow";
  TString r9s = (r9 > 0.94 ? "highR9" : "lowR9");
  
  if(fabs(eta) > 1.0)  region = "EBHigh";
  if(fabs(eta) > 1.48) region = "EELow";
  if(fabs(eta) > 2.0)  region = "EEHigh";

  return scaleSys[ Form("%s_%s",region.Data(),r9s.Data()) ]; 
}

void DataFitter::fixNorm(TString catName, float norm) {
  fixScales = true;
  normMap[catName] = norm;
}

void DataFitter::SetTriggerPath(TString triggerFilePath) {
  if(triggerFilePath!="") {
    std::cout << "adding trigger path: " << triggerFilePath << std::endl;
    fChain->AddFriend("SusyHggTriggerTree",triggerFilePath);
    hasTrigger=true;
    fChain->SetBranchAddress("passTrigger",&triggerBit);
  }
}

void DataFitter::SetNoisePath(TString noiseFilePath) {
  if(noiseFilePath!="") {
    std::cout << "adding noise path: " << noiseFilePath << std::endl;
    fChain->AddFriend("SusyHggNoiseTree",noiseFilePath);
    hasNoise=true;
    fChain->SetBranchAddress("passNoise",&noiseBit);
  }
}
