#include "MCBackgroundFitterNew.hh"

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

void MCBackgroundFitter::buildSidebandHistograms() {
  for(auto cat: catNames) {
    SidebandRegionHistograms[cat] = new TH2F("data_"+cat+"_SidebandRegion","",nXbins-1,xBins,nYbins-1,yBins);
    SidebandRegionHistogramsFineBin[cat] = new TH2F("data_"+cat+"_SidebandRegion_FineBin","",250,0,2500,100,0,1);

    SidebandRegionHistograms[cat+"_fit_Up"] = new TH2F("data_"+cat+"_SidebandRegion_fitUp","",nXbins-1,xBins,nYbins-1,yBins);

    SidebandRegionHistograms[cat+"_fit_Down"] = new TH2F("data_"+cat+"_SidebandRegion_fitDown","",nXbins-1,xBins,nYbins-1,yBins);

    SidebandRegionHistograms[cat+"_bkgShape_Up"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeUp","",nXbins-1,xBins,nYbins-1,yBins);

    SidebandRegionHistograms[cat+"_bkgShape_Down"] = new TH2F("data_"+cat+"_SidebandRegion_bkgShapeDown","",nXbins-1,xBins,nYbins-1,yBins);

    weightMap[cat]=0;

    TTree massTree("massTree","");
    float m;
    massTree.Branch("mgg",&m);
    Long64_t iEntry=-1;
    while( fChain->GetEntry(++iEntry) ) {
      if(!passBasicSelection()) continue;
      TLorentzVector pho1;
      TLorentzVector pho2;
  
      pho1.SetPtEtaPhiM(pho1_pt,pho1_eta,pho1_phi,0);
      pho2.SetPtEtaPhiM(pho2_pt,pho2_eta,pho2_phi,0);
      
      float se1=pho1_sigEoE;
      float se2=pho2_sigEoE;
      
      float btag = highest_csv;

      if( getCategory(pho1,pho2,se1,se2,btag,mbb_NearH,mbb_NearZ) == cat ) {
	if(mgg > 125 - nSigEffSignalRegion*sigmaEffectives[cat] && mgg < 126 + nSigEffSignalRegion*sigmaEffectives[cat] && MR>=150) {
	  weightMap[cat]++;
	}
	m=mgg;
	massTree.Fill();
      }
    }
    weightMap[cat] = normMap[cat]/weightMap[cat];

    std::cout << cat << std::endl;
    RooWorkspace *ws = new RooWorkspace(cat+"_mgg_workspace","");
    scales[cat] = doFitGetScale(&massTree,nSigEffSignalRegion*sigmaEffectives[cat],ws,(cat=="HighPt" ? false:true) );
    mggFitWorkspaces.push_back(ws);
  }
}

void MCBackgroundFitter::Run() {
  buildSidebandHistograms();
  Long64_t iEntry=-1;
  while(fChain->GetEntry(++iEntry)) {
    if(!passBasicSelection()) continue;
    TLorentzVector pho1;
    TLorentzVector pho2;
  
    pho1.SetPtEtaPhiM(pho1_pt,pho1_eta,pho1_phi,0);
    pho2.SetPtEtaPhiM(pho2_pt,pho2_eta,pho2_phi,0);
      
    float se1=pho1_sigEoE;
    float se2=pho2_sigEoE;
      
    float btag = highest_csv;
    TString cat = getCategory(pho1,pho2,se1,se2,btag,mbb_NearH,mbb_NearZ);
    weight=weightMap[cat];

    processEntry();
    processEntrySidebands();
  }
  std::cout << "nProcessed:  " << iEntry << std::endl;
  for(auto cat: catNames) {
    std::cout << cat << "  " << SidebandRegionHistograms[cat]->Integral() << std::endl;

    SidebandRegionHistograms[cat+"_bkgShape_Up"]->Scale( SidebandRegionHistograms[cat]->Integral()/SidebandRegionHistograms[cat+"_bkgShape_Up"]->Integral() );
    SidebandRegionHistograms[cat+"_bkgShape_Down"]->Scale( SidebandRegionHistograms[cat]->Integral()/SidebandRegionHistograms[cat+"_bkgShape_Down"]->Integral() );
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
}

void MCBackgroundFitter::fixNorm(TString catName, float norm) {
  normMap[catName] = norm;
}
