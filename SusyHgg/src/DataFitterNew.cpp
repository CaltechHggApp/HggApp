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

#include <cmath>

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))

RealVar DataFitter::doFitGetScale(TTree* data,float width,RooWorkspace* ws,bool dExp) {

  assert(data->GetEntries() != 0);

  RooRealVar mgg("mgg","",minMgg,maxMgg);

  mgg.setRange("fit",103,160);
  mgg.setRange("low",103,120);
  mgg.setRange("high",131,160);

  RooFormulaVar mggsq("mggsq","","@0*@0",mgg);

  RooRealVar a1("a1","",0.6,-1.,1.);
  RooRealVar a2("a2","",0.4/150,-1/150.,1/150.);
  RooRealVar f("f","",0.5,0,1);


  RooFormulaVar a1sq("a1sq","","-1*@0*@0",a1);
  RooFormulaVar a2sq("a2sq","","-1*@0*@0",a2);
  RooFormulaVar ftanh("ftanh","","0.5*(tanh(@0)+1)",f);

  RooRealVar NBkg1("Nbkg1","",10,-1e5,1e5);
  RooRealVar NBkg2("Nbkg2","",1,-1e5,1e5);

  RooFormulaVar NBkg1Sq("Nbkg1Sq","","@0*@0",NBkg1);
  RooFormulaVar NBkg2Sq("Nbkg2Sq","","@0*@0",NBkg2);

  RooExponential e1("e1","",mgg,a1sq);
  RooExponential e2("e2","",mgg,a2sq);

  RooAddPdf pdf("pdf","",RooArgSet(e1,e2),RooArgSet(NBkg1Sq,NBkg2Sq));

  RooDataSet rdata("data","",data,mgg);
  rdata.Print();
  //pdf_FULL->fitTo(rdata,RooFit::Strategy(0),RooFit::Extended(kTRUE));
  //pdf_FULL->fitTo(rdata,RooFit::Strategy(2),RooFit::Extended(kTRUE),RooFit::Save(kTRUE));

  //pdf.fitTo(rdata,RooFit::Strategy(0),RooFit::Extended(kTRUE));
  //pdf.fitTo(rdata,RooFit::Strategy(2),RooFit::Extended(kTRUE));
  
  //full fit
  pdf.fitTo(rdata,RooFit::Strategy(0),RooFit::Extended(kTRUE));
  RooFitResult* res = pdf.fitTo(rdata,RooFit::Strategy(2),RooFit::Extended(kTRUE),RooFit::Save(kTRUE));
  
  //pdf->fitTo(rdata,RooFit::Strategy(0),RooFit::Extended(kTRUE),RooFit::Minimizer("Minuit2"));
  //RooFitResult* res = pdf->fitTo(rdata,RooFit::Strategy(2),RooFit::Extended(kTRUE),RooFit::Save(kTRUE),RooFit::Minimizer("Minuit2"));
  
  //sideband fit mgg in (103-120, 131-160) GeV
  //pdf.fitTo(rdata,RooFit::Strategy(0),RooFit::Extended(kTRUE),RooFit::Range("low,high"));
  //RooFitResult* res = pdf.fitTo(rdata,RooFit::Strategy(2),RooFit::Extended(kTRUE),RooFit::Save(kTRUE),RooFit::Range("low,high"));
  
  
  std::cout << "======================================" << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "[INFO] mgg width: " << width << " GeV " << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "======================================" << std::endl;

  mgg.setRange("sig",125-width,126+width);

  RooAbsReal *sig_int = pdf.createIntegral(mgg,RooFit::NormSet(mgg),RooFit::Range("sig"));
  sig_int->SetName("signal_integral");

  RooAbsReal *total_int = pdf.createIntegral(mgg);
  total_int->SetName("total_integral");
  
  RooAbsReal *sig_unorm_int = pdf.createIntegral( mgg, RooFit::Range("sig") );
  sig_unorm_int->SetName("sig_unorm_int");

  /*
  std::cout << "======================================" << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "signal_integral: " << sig_int->getVal() << std::endl;
  std::cout << "total_integral: "  << total_int->getVal() << std::endl;
  std::cout << "signal_unorm_integral: "  << sig_unorm_int->getVal() << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "======================================" << std::endl;
  */
  float N_sideband = rdata.sumEntries(Form("(mgg>%0.2f && mgg <120) || (mgg>131 && mgg<%0.2f)",minMgg,maxMgg));  

  RealVar scale;

  /*
  std::cout << "======================================" << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "fitted norm: " << NBkg1Sq.getVal()+NBkg2Sq.getVal() << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "======================================" << std::endl;
  */
  scale.val = (NBkg1Sq.getVal()+NBkg2Sq.getVal())*sig_int->getVal()/N_sideband;
  scale.error = scale.val*sqrt( 1/N_sideband + pow(sig_int->getPropagatedError(*res)/sig_int->getVal(),2) );

  if(isnan(scale.val) || isnan(scale.error)) {
    std::cerr << "scale val: " << scale.val << "   error: " << scale.error << std::endl;
    assert(false);
  }
  RooRealVar scaleFactor("scaleFactor","",scale.val);
  scaleFactor.setError(scale.error);

  if(ws) {
    ws->import(pdf);
    ws->import(rdata);
    ws->import(*res);
    ws->import(*sig_int);
    ws->import(scaleFactor);
  }else{
    delete res;
    delete sig_int;
  }

  std::cout << "=================================" << std::endl;
  std::cout << "=================================" << std::endl;
  std::cout << "[INFO] scaleFactor: " << scale.val << std::endl;
  std::cout << "=================================" << std::endl;
  std::cout << "=================================" << std::endl;
  
  return scale;
}

void DataFitter::buildSidebandHistograms() {
  for(auto cat: catNames) {
    if(useHT) {
      SidebandRegionHistograms[cat] = new TH2F("data_"+cat+"_SidebandRegion","",47,150,2500,12,0,300);
      SidebandRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SidebandRegion_FineBin","",500,0,2500,100,0,300);
      SidebandRegionHistograms[cat+"_statistics_Up"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsUp","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_statistics_Down"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsDown","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_statistics"] = new TH2F("data_"+cat+"_SidebandRegion_statistics","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_statistics_High"] = new TH2F("data_"+cat+"_SidebandRegion_statistics_High","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_statistics_Low"] = new TH2F("data_"+cat+"_SidebandRegion_statistics_Low","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_fit_Up"] = new TH2F("data_"+cat+"_SidebandRegion_fitUp","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_fit_Down"] = new TH2F("data_"+cat+"_SidebandRegion_fitDown","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_sfSyst_Up"] = new TH2F("data_"+cat+"_SidebandRegion_sfSystUp","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_sfSyst_Down"] = new TH2F("data_"+cat+"_SidebandRegion_sfSystDown","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_bkgShapeReal_Up"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeRealUp","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_bkgShapeReal_Down"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeRealDown","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_bkgShape_Up"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeUp","",47,150,2500,12,0,300);
      SidebandRegionHistograms[cat+"_bkgShape_Down"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeDown","",47,150,2500,12,0,300);

      SidebandRegions3D[cat] = new TH3F("data_"+cat+"_SidebandRegion_3D","",47,150,2500,12,0,300,200,100,200);
    }else{
      SidebandRegionHistograms[cat] = new TH2F("data_"+cat+"_SidebandRegion","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SidebandRegion_FineBin","",250,0,2500,100,0,1);
      SidebandRegionHistograms[cat+"_statistics_Up"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_statistics_Down"] = new TH2F("data_"+cat+"_SidebandRegion_statisticsDown","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_statistics"] = new TH2F("data_"+cat+"_SidebandRegion_statistics","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_statistics_High"] = new TH2F("data_"+cat+"_SidebandRegion_statistics_High","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_statistics_Low"] = new TH2F("data_"+cat+"_SidebandRegion_statistics_Low","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_fit_Up"] = new TH2F("data_"+cat+"_SidebandRegion_fitUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_fit_Down"] = new TH2F("data_"+cat+"_SidebandRegion_fitDown","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_sfSyst_Up"] = new TH2F("data_"+cat+"_SidebandRegion_sfSystUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_sfSyst_Down"] = new TH2F("data_"+cat+"_SidebandRegion_sfSystDown","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_bkgShape_Up"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_bkgShape_Down"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeDown","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_bkgShapeReal_Up"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeRealUp","",nXbins-1,xBins,nYbins-1,yBins);
      SidebandRegionHistograms[cat+"_bkgShapeReal_Down"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeRealDown","",nXbins-1,xBins,nYbins-1,yBins);

      SidebandRegions3D[cat] = new TH3F("data_"+cat+"_SidebandRegion_3D","",nXbins-1,xBins,nYbins-1,yBins,nZbins-1,zBins);
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
    //scales[cat] = doFitGetScale(&massTree,nSigEffSignalRegion*sigmaEffectives[cat],ws,(cat=="HighPt" ? false:true) );
    scales[cat] = doFitGetScale(&massTree,nSigEffSignalRegion*sigmaEffectives[cat],ws,true);

    mggFitWorkspaces.push_back(ws);
  }
}

void DataFitter::Run() {
  std::cout << "\n\n\n\n-------------------\n";
  std::cout << "RUNNING ON DATA" << std::endl;
  std::cout << "-------------------\n\n\n" << std::endl;

  for( auto cat: catNames) {
    nSideband[cat]=0;
    nSidebandLow[cat]=0;
    nSidebandHigh[cat]=0;    
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
    assert(!isnan(weight));
    processEntry(false);
    assert(!isnan(weight));
    processEntrySidebands();
  }
  std::cout << "nProcessed:  " << iEntry << std::endl;
  for(auto cat: catNames) {
    std::cout << cat << "  " << SidebandRegionHistograms[cat]->Integral() << std::endl;

    SidebandRegionHistograms[cat+"_bkgShapeReal_Up"]->Scale( nSideband[cat]/nSidebandHigh[cat] ); 
    SidebandRegionHistograms[cat+"_bkgShapeReal_Down"]->Scale( nSideband[cat]/nSidebandLow[cat] ); 
  }
  
  for(auto cat: catNames) {
    assert(nSideband[cat]!=0);
    //float scale = nSideband[cat]/SidebandRegionHistograms[cat]->Integral(); //scale up to full statistics
    for(int xBin=0; xBin<SidebandRegionHistograms[cat]->GetNbinsX()+1; xBin++) {
      for(int yBin=0; yBin<SidebandRegionHistograms[cat]->GetNbinsY()+1; yBin++) {
	float content = SidebandRegionHistograms[cat]->GetBinContent(xBin,yBin);
	float error = TMath::Sqrt(content*scales[cat].val);
	assert(!isnan(error));
	SidebandRegionHistograms[ cat+"_statistics_Up" ]->SetBinContent(xBin,yBin,content+error);
	SidebandRegionHistograms[ cat+"_statistics_Down" ]->SetBinContent(xBin,yBin,content-error);      
      }
    }
  }

  for(auto cat: catNames) {
    for(int xBin=1; xBin<SidebandRegionHistograms[cat]->GetNbinsX()+1; xBin++) {
      for(int yBin=1; yBin<SidebandRegionHistograms[cat]->GetNbinsY()+1; yBin++) {
	float statUp = SidebandRegionHistograms[ cat+"_statistics_Up" ]->GetBinContent(xBin,yBin);
	float statDown = SidebandRegionHistograms[ cat+"_statistics_Down" ]->GetBinContent(xBin,yBin);
	SidebandRegionHistograms[cat+"_bkgShape_Up"]->SetBinContent(xBin,yBin,MAX(SidebandRegionHistograms[cat+"_bkgShapeReal_Up"]->GetBinContent(xBin,yBin),statUp));
	SidebandRegionHistograms[cat+"_bkgShape_Down"]->SetBinContent(xBin,yBin,MIN(SidebandRegionHistograms[cat+"_bkgShapeReal_Down"]->GetBinContent(xBin,yBin),statDown));
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
  for(auto hist: SidebandRegions3D) hist.second->Write();

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
    float thisRsq = Rsq;//Rsq points to t1Rsq, see L493 @ SusyHggTreeBase.h

    if(useHT) {
      thisMR = HT;
      thisRsq = MET;
    }

    if(fabs(weight-1)>0.0001) {
      std::cerr << "DATA has weight: " << weight << std::endl;
      assert(weight==1);
    }
    
    /*
    std::cout << "[INFO]: scaleFactor: " << scales[cat].val << std::endl;
    std::cout << "[INFO]: weight: " << weight << std::endl;
    std::cout << "[INFO]: cat: " << cat << std::endl;
    std::cout << "[INFO]: minMgg: " << minMgg << " maxMgg: " << maxMgg << std::endl;
    */

    assert(!isnan(weight));
    if((mgg>minMgg && mgg<120) || (mgg>131 && mgg<maxMgg)){
      SidebandRegionHistograms[cat]->Fill(thisMR,thisRsq,scales[cat].val*weight);
      SidebandRegionHistograms[cat+"_statistics"]->Fill(thisMR,thisRsq,weight*weight);
      SidebandRegions3D[cat]->Fill(thisMR,thisRsq,mgg,scales[cat].val*weight);
      SidebandRegionHistograms[cat+"_fit_Up"]->Fill(thisMR,thisRsq,(scales[cat].val+scales[cat].error)*weight);
      SidebandRegionHistograms[cat+"_fit_Down"]->Fill(thisMR,thisRsq,(scales[cat].val-scales[cat].error)*weight);
      SidebandRegionHistograms[cat+"_sfSyst_Up"]->Fill(thisMR,thisRsq,(scales[cat].val*(1+sfSyst[cat]))*weight);
      SidebandRegionHistograms[cat+"_sfSyst_Down"]->Fill(thisMR,thisRsq,(scales[cat].val*(1-sfSyst[cat]))*weight);
      SidebandRegionHistogramsFineBin[cat]->Fill(thisMR,thisRsq,scales[cat].val*weight);
      nSideband[cat]++;
    }
    if(mgg > 131 && mgg < maxMgg) {
      SidebandRegionHistograms[cat+"_bkgShapeReal_Up"]->Fill(thisMR,thisRsq,scales[cat].val*weight);
      nSidebandHigh[cat]++;
      SidebandRegionHistograms[cat+"_statistics_High"]->Fill(thisMR,thisRsq,weight*weight);
    }
    if(mgg > minMgg && mgg < 120) {
      SidebandRegionHistograms[cat+"_bkgShapeReal_Down"]->Fill(thisMR,thisRsq,scales[cat].val*weight);
      nSidebandLow[cat]++;
      SidebandRegionHistograms[cat+"_statistics_Low"]->Fill(thisMR,thisRsq,weight*weight);
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
