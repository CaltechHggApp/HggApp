
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"

#include <iostream>

std::pair<float,float> getFitRange(TH1D& hist) {
  int maxBin = hist.GetMaximumBin();
  float targetVal = hist.GetBinContent( maxBin ) /2.;


  int nBins = hist.GetNbinsX();
  
  std::pair<float,float> results = {hist.GetBinLowEdge(1),hist.GetBinLowEdge(nBins+1)};

  int diff=0;
  while( maxBin-(++diff) >0) {
    if( hist.GetBinContent(maxBin-diff) < targetVal ){
      results.first = hist.GetBinCenter(maxBin-diff);
      break;
    }
  }

  diff=0;
  while( maxBin+(++diff) <=nBins) {
    if( hist.GetBinContent(maxBin+diff) < targetVal ){
      results.second = hist.GetBinCenter(maxBin+diff);
      break;
    }
  }
  
  
  return results;
}

void MeasureSmearing(TTree *Full, TTree *Fast,TString tag="") {
  //TChain *Full = loadList("/home/amott/HggApp/Hgg_53X/Caltech/Reduced/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_STAR53_V7A-v1.list","HggReduce");
  //TChain *Fast = loadList("/home/amott/HggApp/Hgg_53X/Caltech/Reduced/SMS-TChiHH_2b2g_2J_mChargino-130to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola__Summer12-START53_V19_FSIM-v1.list","HggReduce");
  

  //TString var = "(Photons.energy-Photons.genMatch.energy)/Photons.genMatch.energy";
  TString var_pho1 = "(pho1_pt*cosh(pho1_eta)-pho1_energyGen)/pho1_energyGen";
  TString var_pho2 = "(pho2_pt*cosh(pho2_eta)-pho2_energyGen)/pho2_energyGen";

  const int nCats=8;

//   TString cats[nCats] = {
//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)<1 && Photons.SC.r9>0.94",
//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1 && abs(Photons.SC.eta)<1.4442 && Photons.SC.r9>0.94",
//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1.56 && abs(Photons.SC.eta)<2 && Photons.SC.r9>0.94",
//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1 && abs(Photons.SC.eta)<2.5 && Photons.SC.r9>0.94",

//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)<1 && Photons.SC.r9<0.94",
//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1 && abs(Photons.SC.eta)<1.4442 && Photons.SC.r9<0.94",
//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1.56 && abs(Photons.SC.eta)<2 && Photons.SC.r9<0.94",
//     "Photons.genMatch.index!=-1 && Photons.energy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1 && abs(Photons.SC.eta)<2.5 && Photons.SC.r9<0.94",
//   };

  TString cats_pho1[nCats] = {
    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=0.00 && abs(pho1_eta)<1.000 && pho1_r9 > 0.94",
    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.00 && abs(pho1_eta)<1.442 && pho1_r9 > 0.94",
    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.56 && abs(pho1_eta)<2.000 && pho1_r9 > 0.94",
    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=2.00 && abs(pho1_eta)<2.500 && pho1_r9 > 0.94",

    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=0.00 && abs(pho1_eta)<1.000 && pho1_r9 <= 0.94",
    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.00 && abs(pho1_eta)<1.442 && pho1_r9 <= 0.94",
    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.56 && abs(pho1_eta)<2.000 && pho1_r9 <= 0.94",
    "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=2.00 && abs(pho1_eta)<2.500 && pho1_r9 <= 0.94",
  };

  TString cats_pho2[nCats] = {
    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=0.00 && abs(pho2_eta)<1.000 && pho2_r9 > 0.94",
    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.00 && abs(pho2_eta)<1.442 && pho2_r9 > 0.94",
    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.56 && abs(pho2_eta)<2.000 && pho2_r9 > 0.94",
    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=2.00 && abs(pho2_eta)<2.500 && pho2_r9 > 0.94",

    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=0.00 && abs(pho2_eta)<1.000 && pho2_r9 <= 0.94",
    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.00 && abs(pho2_eta)<1.442 && pho2_r9 <= 0.94",
    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.56 && abs(pho2_eta)<2.000 && pho2_r9 <= 0.94",
    "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=2.00 && abs(pho2_eta)<2.500 && pho2_r9 <= 0.94",
  };

  TString catNames[nCats] = {
    "EBLow_Gold",
    "EBHigh_Gold",
    "EELow_Gold",
    "EEHigh_Gold",
    "EBLow_Bad",
    "EBHigh_Bad",
    "EELow_Bad",
    "EEHigh_Bad",
  };

  TCanvas cv;
  TH1D pt_full("pt_full","",80,0,400);
  pt_full.SetXTitle("p_{T} [GeV]");
  pt_full.SetYTitle("Events / 5 GeV");

  //Full->Project("pt_full","Photons.energy/cosh(Photons.eta)","Photons.genMatch.index != -1");

  TH1D pt_fast("pt_fast","",80,0,400);
  pt_fast.SetXTitle("p_{T} [GeV]");
  pt_fast.SetYTitle("Events / 5 GeV");

  //Fast->Project("pt_fast","Photons.energy/cosh(Photons.eta)","Photons.genMatch.index != -1");

  TH1D* ratio = (TH1D*)pt_fast.Clone("pt_ratio");
  ratio->Divide(&pt_full);



  TH1D h("h","",50,-0.2,0.2);
  TH1D h_p2("h_p2","",50,-0.2,0.2);
  h.SetXTitle("(E_{reco}-E_{MC})/E_{MC}");

  for(int i=0;i<nCats;i++) {    
    std::cout << catNames[i] << std::endl;
    Full->Project("h",var_pho1,cats_pho1[i]);
    Full->Project("h_p2",var_pho2,cats_pho2[i]);
    h.Add(&h_p2);
    std::pair<float,float> r_full = getFitRange(h);
    std::cout << "FullSIM fit range: " << r_full.first << " - " << r_full.second << std::endl;
    h.Fit("gaus","","",r_full.first,r_full.second);
    h.Draw();
    cv.SaveAs("ValidationPlots/Smearing_FullSIM_"+catNames[i]+tag+".png");

    Fast->Project("h",var_pho1,cats_pho1[i]+" && m23 < 175.");
    Fast->Project("h_p2",var_pho2,cats_pho2[i]+" && m23 < 175.");
    h.Add(&h_p2);
    std::pair<float,float> r_fast = getFitRange(h);
    std::cout << "FastSIM fit range: " << r_fast.first << " - " << r_fast.second << std::endl;
    h.Fit("gaus","","",r_fast.first,r_fast.second);
    h.Draw();
    cv.SaveAs("ValidationPlots/Smearing_FastSIM_"+catNames[i]+tag+".png");    
  }
}


void MeasureSmearing(TString fullFileName, TString fastFileName,TString treeName="SusyHggTree",TString tag="") {
  TFile fullFile(fullFileName);
  TTree* full = (TTree*)fullFile.Get(treeName);
  TFile fastFile(fastFileName);
  TTree* fast = (TTree*)fastFile.Get(treeName);

  MeasureSmearing(full,fast,tag);

}
