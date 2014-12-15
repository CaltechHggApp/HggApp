
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"

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

//Crystal ball function for signal, parameters are 0:alpha,1:n,2:mean,3:sigma,4:normalization;

Double_t CrystalBall(Double_t *x,Double_t *par) {

  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);

  if (t >= -absAlpha) {
    return par[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha; 

    return par[4]*(a/TMath::Power(b - t, par[1]));
  }
}

//Double Crystal ball function for signal, parameters are 0:alpha1,1:n1,2:alpha2,3:n2,4:mean,5:sigma,6:normalization;
// alpha1 and alpha2 will be forced positive

Double_t DoubleCrystalBall(Double_t *x,Double_t *par) {

  Double_t t = (x[0]-par[4])/par[5];

  Double_t absAlpha1 = fabs((Double_t)par[0]);
  Double_t absAlpha2 = fabs((Double_t)par[2]);

  if (t >= -absAlpha1 && t <=absAlpha2) {
    return par[6]*exp(-0.5*t*t);
  }
  else if(t<-absAlpha1){
    Double_t a =  TMath::Power(par[1]/absAlpha1,par[1])*exp(-0.5*absAlpha1*absAlpha1);
    Double_t b= par[1]/absAlpha1 - absAlpha1; 

    return par[6]*(a/TMath::Power(b - t, par[1]));
  }else{
    Double_t a =  TMath::Power(par[3]/absAlpha2,par[3])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t b= par[3]/absAlpha2 - absAlpha2; 

    return par[6]*(a/TMath::Power(b + t, par[3]));    
  }
}

void MeasureSmearing(TTree *Full, TTree *Fast,TString tag="",bool useElectrons=false) {
  //TChain *Full = loadList("/home/amott/HggApp/Hgg_53X/Caltech/Reduced/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_STAR53_V7A-v1.list","HggReduce");
  //TChain *Fast = loadList("/home/amott/HggApp/Hgg_53X/Caltech/Reduced/SMS-TChiHH_2b2g_2J_mChargino-130to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola__Summer12-START53_V19_FSIM-v1.list","HggReduce");
  
  Full->SetBranchStatus("*",0);
  if(useElectrons) Full->SetBranchStatus("Electrons.*",1);
  else Full->SetBranchStatus("Photons.*",1);

  Fast->SetBranchStatus("*",0);
  if(useElectrons) Fast->SetBranchStatus("Electrons.*",1);
  else Fast->SetBranchStatus("Photons.*",1);


  TString var = "Photons.correctedEnergy/Photons.genMatch.energy-1";
  if(useElectrons) var = "Electrons.correctedEnergy/Electrons.genMatch.energy-1";

  TString epVar = "Electrons.EOverP";
//   TString var_pho1 = "(pho1_pt*cosh(pho1_eta)-pho1_energyGen)/pho1_energyGen";
//   TString var_pho2 = "(pho2_pt*cosh(pho2_eta)-pho2_energyGen)/pho2_energyGen";

  TString FSIM_SCALE = "( (abs(Photons.SC.eta)<1.56)*1 + (abs(Photons.SC.eta)>=1.56 && abs(Photons.SC.eta)<2)*( (Photons.SC.r9>0.94)*1.0276 + (Photons.SC.r9<=0.94)*1.0195) + (abs(Photons.SC.eta)>=2. && abs(Photons.SC.eta)<2.5)*( (Photons.SC.r9>0.94)*1.0216 + (Photons.SC.r9<=0.94)*1.0330) )*";

  const int nCats=8;

   TString cats[nCats] = {
     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)<1 && Photons.SC.r9>0.94",
     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1 && abs(Photons.SC.eta)<1.4442 && Photons.SC.r9>0.94",
     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1.56 && abs(Photons.SC.eta)<2 && Photons.SC.r9>0.94",
     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=2 && abs(Photons.SC.eta)<2.5 && Photons.SC.r9>0.94",

     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)<1 && Photons.SC.r9<0.94",
     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1 && abs(Photons.SC.eta)<1.4442 && Photons.SC.r9<0.94",
     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=1.56 && abs(Photons.SC.eta)<2 && Photons.SC.r9<0.94",
     "Photons.genMatch.index!=-1 && Photons.correctedEnergy/cosh(Photons.eta)>25 && abs(Photons.SC.eta)>=2 && abs(Photons.SC.eta)<2.5 && Photons.SC.r9<0.94",
   };

   if(useElectrons) {
     FSIM_SCALE.ReplaceAll("Photons","Electrons");
     for(int i=0;i<nCats;i++) {
       cats[i].ReplaceAll("Photons","Electrons");
       std::cout << cats[i] << std::endl;
     }
   }

//   TString cats_pho1[nCats] = {
//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=0.00 && abs(pho1_eta)<1.000 && pho1_r9 > 0.94",
//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.00 && abs(pho1_eta)<1.442 && pho1_r9 > 0.94",
//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.56 && abs(pho1_eta)<2.000 && pho1_r9 > 0.94",
//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=2.00 && abs(pho1_eta)<2.500 && pho1_r9 > 0.94",

//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=0.00 && abs(pho1_eta)<1.000 && pho1_r9 <= 0.94",
//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.00 && abs(pho1_eta)<1.442 && pho1_r9 <= 0.94",
//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=1.56 && abs(pho1_eta)<2.000 && pho1_r9 <= 0.94",
//     "pho1_genMatch==1 && pho1_pt>25 && abs(pho1_eta)>=2.00 && abs(pho1_eta)<2.500 && pho1_r9 <= 0.94",
//   };

//   TString cats_pho2[nCats] = {
//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=0.00 && abs(pho2_eta)<1.000 && pho2_r9 > 0.94",
//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.00 && abs(pho2_eta)<1.442 && pho2_r9 > 0.94",
//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.56 && abs(pho2_eta)<2.000 && pho2_r9 > 0.94",
//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=2.00 && abs(pho2_eta)<2.500 && pho2_r9 > 0.94",

//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=0.00 && abs(pho2_eta)<1.000 && pho2_r9 <= 0.94",
//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.00 && abs(pho2_eta)<1.442 && pho2_r9 <= 0.94",
//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=1.56 && abs(pho2_eta)<2.000 && pho2_r9 <= 0.94",
//     "pho2_genMatch==1 && pho2_pt>25 && abs(pho2_eta)>=2.00 && abs(pho2_eta)<2.500 && pho2_r9 <= 0.94",
//   };

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



  TH1D h("h","",100,-0.2,0.2);
  TH1D h_p2("h_p2","",100,-0.2,0.2);
  h.SetXTitle("(E_{reco}-E_{MC})/E_{MC}");

  TH1D ep("ep","",200,0.6,1.4);
  ep.SetXTitle("E/p");

  for(int i=2;i<nCats;i++) {    
    std::cout << catNames[i] << std::endl;
    Full->Project("h",var,cats[i]);
//     Full->Project("h",var_pho1,cats_pho1[i]);
//     Full->Project("h_p2",var_pho2,cats_pho2[i]);
    //h.Add(&h_p2);
    std::pair<float,float> r_full = getFitRange(h);
    std::cout << "FullSIM fit range: " << r_full.first << " - " << r_full.second << std::endl;
    TF1 full_fit("full_fit",DoubleCrystalBall,-0.2,0.2,7);
    full_fit.SetParNames("#alpha_{1}","n_{1}","#alpha_{2}","n_{2}","Mean","#sigma","N");
    full_fit.SetParLimits(0,0.2,20);
    full_fit.SetParLimits(1,0.00001,100000);
    full_fit.SetParLimits(2,0.2,20);
    full_fit.SetParLimits(3,0.00001,100000);
    full_fit.SetParLimits(4,-0.1,0.1);
    full_fit.SetParLimits(5,0.001,1);
    full_fit.SetParameter(5,0.01);
    full_fit.SetParLimits(6,100,1e9);

    full_fit.SetParameter(0,1);
    full_fit.SetParameter(1,4);
    full_fit.SetParameter(2,1);
    full_fit.SetParameter(3,2);
    full_fit.SetParameter(4,0);
    full_fit.SetParameter(6,10000);
    //full_fit.FixParameter(6,h.Integral());

    h.Fit(&full_fit,"VEM");
    //h.Fit("gaus","","",r_full.first,r_full.second);
    h.Draw();
    cv.SaveAs("ValidationPlots/Smearing_FullSIM_"+catNames[i]+tag+".png");
    
    if(useElectrons && false) {
      Full->Project("ep",epVar,cats[i]);
      std::pair<float,float> r_ep_full = getFitRange(ep);
      std::cout << "FullSIM fit range: " << r_ep_full.first << " - " << r_ep_full.second << std::endl;      
      TF1 full_fit_ep("full_fit_ep",DoubleCrystalBall,0.6,1.4,7);
      full_fit_ep.SetParNames("#alpha_{1}","n_{1}","#alpha_{2}","n_{2}","Mean","#sigma","N");
      full_fit_ep.SetParLimits(0,0.2,20);
      full_fit_ep.SetParLimits(1,0.00001,100000);
      full_fit_ep.SetParLimits(2,0.2,20);
      full_fit_ep.SetParLimits(3,0.00001,100000);
      full_fit_ep.SetParLimits(4,0.9,1.1);
      full_fit_ep.SetParLimits(5,0.001,1);
      full_fit_ep.SetParameter(5,0.01);
      full_fit_ep.SetParLimits(6,100,1e9);
      full_fit_ep.SetParameter(6,h.Integral());

      full_fit_ep.SetParameter(0,1);
      full_fit_ep.SetParameter(1,4);
      full_fit_ep.SetParameter(2,1);
      full_fit_ep.SetParameter(3,2);
      full_fit_ep.SetParameter(4,1);
      full_fit_ep.SetParameter(6,10000);


      ep.Fit(&full_fit_ep,"VEMI");
      //ep.Fit("gaus","","",r_ep_full.first,r_ep_full.second);
      ep.Draw();
      cv.SaveAs("ValidationPlots/EP_FullSIM_"+catNames[i]+tag+".png");      
    }

    Fast->Project("h",var,cats[i]);
//     Fast->Project("h",var_pho1,cats_pho1[i]+" && m23 < 177.");
//     Fast->Project("h_p2",var_pho2,cats_pho2[i]+" && m23 < 175.");
    //h.Add(&h_p2);
    std::pair<float,float> r_fast = getFitRange(h);
    std::cout << "FastSIM fit range: " << r_fast.first << " - " << r_fast.second << std::endl;
    TF1 fast_fit("fast_fit",DoubleCrystalBall,-0.2,0.2,7);
    fast_fit.SetParNames("#alpha_{1}","n_{1}","#alpha_{2}","n_{2}","Mean","#sigma","N");
    fast_fit.SetParLimits(0,0.2,20);
    fast_fit.SetParLimits(1,0.00001,100000);
    fast_fit.SetParLimits(2,0.2,20);
    fast_fit.SetParLimits(3,0.00001,100000);
    fast_fit.SetParLimits(4,-0.08,0.08);
    fast_fit.SetParLimits(5,0.001,1);
    fast_fit.SetParameter(5,0.01);
    fast_fit.SetParLimits(6,100,1e9);
    fast_fit.SetParameter(6,h.Integral());

    fast_fit.SetParameter(0,1);
    fast_fit.SetParameter(1,4);
    fast_fit.SetParameter(2,1);
    fast_fit.SetParameter(3,2);
    fast_fit.SetParameter(4,0);
    fast_fit.SetParameter(6,10000);

    h.Fit(&fast_fit,"VEMI");
    //h.Fit("gaus","","",r_fast.first,r_fast.second);
    h.Draw();
    cv.SaveAs("ValidationPlots/Smearing_FastSIM_"+catNames[i]+tag+".png");    

    if(useElectrons) {
      Fast->Project("ep",epVar,cats[i]);
      std::pair<float,float> r_ep_fast = getFitRange(ep);
      std::cout << "FastSIM fit range: " << r_ep_fast.first << " - " << r_ep_fast.second << std::endl;      
      Full->Project("ep",epVar,cats[i]);
      std::pair<float,float> r_ep_full = getFitRange(ep);
      std::cout << "FullSIM fit range: " << r_ep_full.first << " - " << r_ep_full.second << std::endl;      

      TF1 fast_fit_ep("fast_fit_ep",DoubleCrystalBall,0.6,1.4,7);
      fast_fit_ep.SetParNames("#alpha_{1}","n_{1}","#alpha_{2}","n_{2}","Mean","#sigma","N");
      fast_fit_ep.SetParLimits(0,0.2,20);
      fast_fit_ep.SetParLimits(1,0.00001,100000);
      fast_fit_ep.SetParLimits(2,0.2,20);
      fast_fit_ep.SetParLimits(3,0.00001,100000);
      fast_fit_ep.SetParLimits(4,0.9,1.1);
      fast_fit_ep.SetParLimits(5,0.001,1);
      fast_fit_ep.SetParameter(5,0.01);
      fast_fit_ep.SetParLimits(6,100,1e9);
      fast_fit_ep.SetParameter(6,h.Integral());

      fast_fit_ep.SetParameter(0,1);
      fast_fit_ep.SetParameter(1,4);
      fast_fit_ep.SetParameter(2,1);
      fast_fit_ep.SetParameter(3,2);
      fast_fit_ep.SetParameter(4,1);
      fast_fit_ep.SetParameter(6,10000);


      ep.Fit(&fast_fit_ep,"VEMI");
      //ep.Fit("gaus","","",r_ep_fast.first,r_ep_fast.second);
      ep.Draw();
      cv.SaveAs("ValidationPlots/EP_FastSIM_"+catNames[i]+tag+".png");      
    }


    Fast->Project("h",FSIM_SCALE+var,cats[i]);
//     Fast->Project("h",var_pho1,cats_pho1[i]+" && m23 < 177.");
//     Fast->Project("h_p2",var_pho2,cats_pho2[i]+" && m23 < 175.");
    //h.Add(&h_p2);
    std::pair<float,float> r_fast_scale = getFitRange(h);
    std::cout << "Fast_ScaleSIM fit range: " << r_fast_scale.first << " - " << r_fast_scale.second << std::endl;
    TF1 fast_scale_fit("fast_scale_fit",DoubleCrystalBall,-0.2,0.2,7);
    fast_scale_fit.SetParNames("#alpha_{1}","n_{1}","#alpha_{2}","n_{2}","Mean","#sigma","N");
    fast_scale_fit.SetParLimits(0,0.2,20);
    fast_scale_fit.SetParLimits(1,0.00001,100000);
    fast_scale_fit.SetParLimits(2,0.2,20);
    fast_scale_fit.SetParLimits(3,0.00001,100000);
    fast_scale_fit.SetParLimits(4,-0.08,0.08);
    fast_scale_fit.SetParLimits(5,0.001,1);
    fast_scale_fit.SetParameter(5,0.01);
    fast_scale_fit.SetParLimits(6,100,1e9);
    fast_scale_fit.SetParameter(6,h.Integral());

    fast_scale_fit.SetParameter(0,1);
    fast_scale_fit.SetParameter(1,4);
    fast_scale_fit.SetParameter(2,1);
    fast_scale_fit.SetParameter(3,2);
    fast_scale_fit.SetParameter(4,0);
    fast_scale_fit.SetParameter(6,10000);

    h.Fit(&fast_scale_fit,"VEMI");
    //h.Fit("gaus","","",r_fast_scale.first,r_fast_scale.second);
    h.Draw();
    cv.SaveAs("ValidationPlots/Smearing_FastSIM_SCALED_"+catNames[i]+tag+".png");    

    if(useElectrons) {
      Fast->Project("ep",epVar,cats[i]);
      std::pair<float,float> r_ep_fast_scale = getFitRange(ep);
      std::cout << "Fast_ScaleSIM fit range: " << r_ep_fast_scale.first << " - " << r_ep_fast_scale.second << std::endl;      
      Full->Project("ep",epVar,cats[i]);
      std::pair<float,float> r_ep_full = getFitRange(ep);
      std::cout << "FullSIM fit range: " << r_ep_full.first << " - " << r_ep_full.second << std::endl;      

      TF1 fast_scale_fit_ep("fast_scale_fit_ep",DoubleCrystalBall,0.6,1.4,7);
      fast_scale_fit_ep.SetParNames("#alpha_{1}","n_{1}","#alpha_{2}","n_{2}","Mean","#sigma","N");
      fast_scale_fit_ep.SetParLimits(0,0.2,20);
      fast_scale_fit_ep.SetParLimits(1,0.00001,100000);
      fast_scale_fit_ep.SetParLimits(2,0.2,20);
      fast_scale_fit_ep.SetParLimits(3,0.00001,100000);
      fast_scale_fit_ep.SetParLimits(4,0.9,1.1);
      fast_scale_fit_ep.SetParLimits(5,0.001,1);
      fast_scale_fit_ep.SetParameter(5,0.01);
      fast_scale_fit_ep.SetParLimits(6,100,1e9);
      fast_scale_fit_ep.SetParameter(6,h.Integral());

      fast_scale_fit_ep.SetParameter(0,1);
      fast_scale_fit_ep.SetParameter(1,4);
      fast_scale_fit_ep.SetParameter(2,1);
      fast_scale_fit_ep.SetParameter(3,2);
      fast_scale_fit_ep.SetParameter(4,1);
      fast_scale_fit_ep.SetParameter(6,10000);

      ep.Fit(&fast_scale_fit_ep);
      //ep.Fit("gaus","","",r_ep_fast_scale.first,r_ep_fast_scale.second);
      ep.Draw();
      cv.SaveAs("ValidationPlots/EP_FastSIM_SCALED_"+catNames[i]+tag+".png");      
    }


  }
}


void MeasureSmearing(TString fullFileName, TString fastFileName,TString treeName="SusyHggTree",TString tag="") {
  TFile fullFile(fullFileName);
  TTree* full = (TTree*)fullFile.Get(treeName);
  TFile fastFile(fastFileName);
  TTree* fast = (TTree*)fastFile.Get(treeName);

  MeasureSmearing(full,fast,tag);

}
