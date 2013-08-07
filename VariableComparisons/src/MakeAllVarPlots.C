#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TTreeFormula.h"
#include "TLatex.h"
#include "TList.h"

#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <array>
#include <assert.h>


#include <cstdlib>

//local includes
#include "include/weightManager.hh"
#include "include/varCorrector.hh"
#include "include/plotManager.hh"


void setstyle();

TCanvas *makeCanvas(std::array<TH1F*,3> data,std::array<TH1F*,3> mc, TString xName,TString label="");

#define __USE_SHERPA  0

int main(int argc,char** argv){
  bool useSherpa=false;
  if(argc>1) {
    if( strcmp(argv[1],"SHERPA")==0 ) useSherpa=true;
  }

  if(useSherpa) std::cout << "SHERPA!!" << std::endl;

  float lumi=6.3;
  weightManager weights;

  plotManager mcPlotter("mc");
  plotManager dataPlotter("data");
  mcPlotter.setUse4Cat();
  dataPlotter.setUse4Cat();

  setstyle();

  std::vector<TString> samples = {
    "GJets_Pt20to40.root",
    "GJets_Pt40.root",
    "QCD_Pt30to40.root",
    "QCD_Pt40.root",
    "DYJetsToLL_M-50.root" };


  if(useSherpa) {
    samples.push_back("DiPhotonJets_sherpa.root");
  }
  else {
    samples.push_back("DiPhotonBox_Pt10to25.root");
    samples.push_back("DiPhotonBox_Pt25to250.root");
    samples.push_back("DiPhotonBox_Pt250.root");
    samples.push_back("DiPhotonJets.root");    
  }

  const int Nsamples = samples.size();

  std::vector<TString> data = {
    "DoublePhoton_Run2012B_13Jul2012.root"
  };

  const int Ndata = data.size();


  std::vector<TString> globalVetos;
  globalVetos.push_back("etaSC > -1.775 && etaSC < -1.764 && phi > 1.2 && phi < 1.6");
  globalVetos.push_back("etaSC > 1.591 && etaSC < 1.595 && phi > -2.06 && phi < -2.045");
  //globalVetos.push_back("etaSC > -0.9 && etaSC < -0.7 && phi > 2.9 && phi < 3.1");
  globalVetos.push_back("etaSC > 1.74 && etaSC < 1.76 && phi > 2.1 && phi < 2.15");
  globalVetos.push_back("etaSC > 1.564 && etaSC < 1.565 && phi > 0.528 && phi < 0.532");

  //output->Draw("phi:etaSC","!(etaSC > -1.8 && etaSC < -1.76 && phi > 1.2 && phi < 1.5) && !(etaSC>1.591 && etaSC<1.595 && phi > -2.06 && phi < -2.045) && !(etaSC > 1.74 && etaSC < 1.76 && phi > 2.1 && phi < 2.15) && !(etaSC > 1.564 && etaSC < 1.565 && phi > 0.528 && phi < 0.532) && !(etaSC > -1.6888 && etaSC < -1.6868 && phi > 1.8 && phi < 1.802) && sieie < 1e-5 && abs(etaSC) > 1.557 && abs(etaSC)<2

  const int Ncat=52;
  TString cats[Ncat] = {
    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==0",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==0",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==0",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==0",
    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==1",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==1",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==1",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 80 && electronMatch==1",

    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==0",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==0",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==0",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==0",
    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==1",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==1",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==1",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 87 && mass < 94 && electronMatch==1",

    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",

    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==0 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==1",

    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0",
    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1",

    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0 && se < 0.012",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0 && se < 0.018",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0 && se < 0.025",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==0 && se < 0.023",
    "abs(etaSC)<1.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1  && se < 0.012",
    "abs(etaSC)>=1.00 && abs(etaSC)<1.44 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1 && se < 0.018",
    "abs(etaSC)>1.57 && abs(etaSC)<2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1  && se < 0.025",
    "abs(etaSC)>=2.00 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 180 && electronMatch==1 && se < 0.023",

    "abs(etaSC)>=2.00 && abs(etaSC) < 2.10 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>=2.10 && abs(etaSC) < 2.20 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>=2.20 && abs(etaSC) < 2.30 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",
    "abs(etaSC)>=2.30 && abs(etaSC) < 2.40 && passPre==1 && passCiC_iso==1 && pt>25 && pt/mass > 1./3. && Trigger==1 && mass > 110 && mass < 150 && electronMatch==0",

  };

  TString catNames[Ncat] = {
    "EBlow_pt25_passCiCIso_pho_Trigger","EBhigh_pt25_passCiCIso_pho_Trigger","EElow_pt25_passCiCIso_pho_Trigger","EEhigh_pt25_passCiCIso_pho_Trigger",
    "EBlow_pt25_passCiCIso_ele_Trigger","EBhigh_pt25_passCiCIso_ele_Trigger","EElow_pt25_passCiCIso_ele_Trigger","EEhigh_pt25_passCiCIso_ele_Trigger",

    "EBlow_pt25_passCiCIso_pho_Z","EBhigh_pt25_passCiCIso_pho_Z","EElow_pt25_passCiCIso_pho_Z","EEhigh_pt25_passCiCIso_pho_Z",
    "EBlow_pt25_passCiCIso_ele_Z","EBhigh_pt25_passCiCIso_ele_Z","EElow_pt25_passCiCIso_ele_Z","EEhigh_pt25_passCiCIso_ele_Z",

    "EBlow_pt25_passCiCIso_pho_110_150","EBhigh_pt25_passCiCIso_pho_110_150","EElow_pt25_passCiCIso_pho_110_150","EEhigh_pt25_passCiCIso_pho_110_150",
    "EBlow_pt25_passCiCIso_ele_110_150","EBhigh_pt25_passCiCIso_ele_110_150","EElow_pt25_passCiCIso_ele_110_150","EEhigh_pt25_passCiCIso_ele_110_150",

    "EBlow_pt25_failCiCIso_pho_110_150","EBhigh_pt25_failCiCIso_pho_110_150","EElow_pt25_failCiCIso_pho_110_150","EEhigh_pt25_failCiCIso_pho_110_150",
    "EBlow_pt25_failCiCIso_ele_110_150","EBhigh_pt25_failCiCIso_ele_110_150","EElow_pt25_failCiCIso_ele_110_150","EEhigh_pt25_failCiCIso_ele_110_150",

    "EBlow_pt25_passCiCIso_pho_110_180","EBhigh_pt25_passCiCIso_pho_110_180","EElow_pt25_passCiCIso_pho_110_180","EEhigh_pt25_passCiCIso_pho_110_180",
    "EBlow_pt25_passCiCIso_ele_110_180","EBhigh_pt25_passCiCIso_ele_110_180","EElow_pt25_passCiCIso_ele_110_180","EEhigh_pt25_passCiCIso_ele_110_180",

    "EBlow_pt25_passCiCIso_se_pho_110_180","EBhigh_pt25_passCiCIso_se_pho_110_180","EElow_pt25_passCiCIso_se_pho_110_180","EEhigh_pt25_passCiCIso_se_pho_110_180",
    "EBlow_pt25_passCiCIso_se_ele_110_180","EBhigh_pt25_passCiCIso_se_ele_110_180","EElow_pt25_passCiCIso_se_ele_110_180","EEhigh_pt25_passCiCIso_se_ele_110_180",

    "EE_2p0_2p1_pt25_passCiCIso_ele_110_150",
    "EE_2p1_2p2_pt25_passCiCIso_ele_110_150",
    "EE_2p2_2p3_pt25_passCiCIso_ele_110_150",
    "EE_2p3_2p4_pt25_passCiCIso_ele_110_150",
  };

  for(int i=0;i<Ncat;i++){
    mcPlotter.addCategory(catNames[i],cats[i]);
    dataPlotter.addCategory(catNames[i],cats[i]);
  }

  mcPlotter.addVetos(&globalVetos);
  dataPlotter.addVetos(&globalVetos);
  const int Nvar=19;
  TString vars[Nvar] = {"se","r9","sieie","sieip","sipip","etaWidth","phiWidth","HE",  "energyBC/rawE","e3x3/energyBC","e5x5/energyBC","eMax/energyBC","e2x5Max/energyBC","nBC","etaBC-etaSC","pt*cosh(etaSC)/rawE","mass","pt","etaSC"};
  int     bins[Nvar] = {50  , 80 , 50 ,   50 ,    50 ,    50 ,       50 ,       50 ,   50 ,            50 ,            50 ,             50 ,           50,                 20,  50,            50,                  80,   100, 50 };
  float   low [Nvar] = {0.0 , 0.2, 0.0,   0.0,    0.0,    0.0,       0.0,       0.0,   0.0,            0.0,            0.0,             0.0,           0.0,                0,   0.0,           0.5,                 80,    0,   -2.5};
  float   high[Nvar] = {0.04, 1.0, 0.04,  0.0004,  0.1,    0.05,      0.2,       0.2,   2.0,            2.0,            2.0,             2.0,           2.0,               20,  0.03,          1.5,                 200,   1000, 2.5};


  for(int i=0;i<Nvar;i++){
    mcPlotter.addVariable  (vars[i],vars[i],bins[i],low[i],high[i]);
    mcPlotter.addVariable  ("corr_"+vars[i],vars[i],bins[i],low[i],high[i],true);
    dataPlotter.addVariable(vars[i],vars[i],bins[i],low[i],high[i]);
    dataPlotter.addVariable("corr_"+vars[i],vars[i],bins[i],low[i],high[i],true);
  }

  for(int iSample=0;iSample<Nsamples;iSample++){
    TFile *f = new TFile("output/"+samples[iSample]);
    TChain *fChain = (TChain*)f->Get("output");
    TString sampleName = samples[iSample];
    sampleName.Remove(sampleName.Last('.'));
    mcPlotter.processChain(fChain,weights.getWeight(sampleName,lumi));
  }

  for(int iData=0;iData<Ndata;iData++){
    TFile *f = new TFile("output/"+data[iData]);
    TChain *fChain = (TChain*)f->Get("output");
    dataPlotter.processChain(fChain,1);
  }


  TString outputFileName = (useSherpa ? "output_histgrams_SHERPA.root" : "output_histograms_pythia.root");
  TFile outputFile(outputFileName,"RECREATE");

  for(int iCat=0;iCat<Ncat;iCat++){
    for(int iVar=0;iVar<Nvar;iVar++){
      std::array<TH1F*,3> mc = mcPlotter.getHistogram(catNames[iCat],vars[iVar]);
      std::array<TH1F*,3> mcCorr = mcPlotter.getHistogram(catNames[iCat],"corr_"+vars[iVar]);
      std::array<TH1F*,3> data = dataPlotter.getHistogram(catNames[iCat],vars[iVar]);
      std::array<TH1F*,3> dataCorr = dataPlotter.getHistogram(catNames[iCat],"corr_"+vars[iVar]);
      TCanvas *cv = makeCanvas(data,mc,vars[iVar],catNames[iCat]);


      TString saveVar = vars[iVar];
      TString folder  = "figs/";
	
      if(useSherpa) {
	saveVar = "SHERPA__"+saveVar;
	folder  = folder+"SHERPA/";
      }
      saveVar.ReplaceAll("/","By");
      saveVar.ReplaceAll("-","_minus_");
      saveVar.ReplaceAll("(","_");
      saveVar.ReplaceAll(")","_");

      
      TString path = folder+saveVar+"_"+catNames[iCat];

      cv->SaveAs(path+"_lin.png");
      ((TPad*)cv->GetPrimitive("plotPad"))->SetLogy();
      cv->SaveAs(path+"_log.png");

      /*
      if(saveVar!="se"){
	//TCanvas *cv2 = makeCanvas(dataCorr,mc);
	TCanvas *cv2 = makeCanvas(data,mcCorr,vars[iVar],catNames[iCat]);
	cv2->SaveAs(Form("%s/%s_corr_%s_lin.png",saveVar.Data(),catNames[iCat].Data()));
	((TPad*)cv2->GetPrimitive("plotPad"))->SetLogy();
	cv2->SaveAs(Form("%s/%s_corr_%s_log.png",saveVar.Data(),catNames[iCat].Data()));

      }else{
	TCanvas *cv2 = makeCanvas(data,mcCorr,vars[iVar],catNames[iCat]);
	cv2->SaveAs(Form("%s/%s_corr_%s_lin.png",saveVar.Data(),catNames[iCat].Data()));
	((TPad*)cv2->GetPrimitive("plotPad"))->SetLogy();
	cv2->SaveAs(Form("%s/%s_corr_%s_log.png",saveVar.Data(),catNames[iCat].Data()));	
      }
      */
      outputFile.cd();
      for(auto h : mc)       h->Write();
      for(auto h : mcCorr)   h->Write();
      for(auto h : data)     h->Write();
      for(auto h : dataCorr) h->Write();
    }
  }  
    outputFile.Close();
    return 0;  
}

void setstyle(){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  
}

TCanvas *makeCanvas(std::array<TH1F*,3> data,std::array<TH1F*,3> mc,TString xName,TString label){
  TH1F *data_total = (TH1F*)data[0]->Clone("data_total");
  data_total->Add(data[1]);
  data_total->Add(data[2]);
  
  float data_norm = data_total->Integral();
  std::cout << "# data events:  " << data_norm << std::endl;
  
  THStack mc_stack("mc_stack","");
  TH1F*   mc_total = (TH1F*)mc[0]->Clone("mc_total");
  mc_total->Add(mc[1]);
  mc_total->Add(mc[2]);
  //compute the total integral for normalization
  double mc_integral=0;
  for(auto h : mc) {
    mc_integral+=h->Integral(); 
  }

  //make the stack
  std::array<Color_t,3> colors = {kBlue,kGreen,kRed};
  assert(colors.size() == mc.size());

  auto h = mc.begin();
  auto c = colors.begin();

  for(; h!=mc.end(); h++,c++) {
    //scale MC to data
    (*h)->Scale(data_norm/mc_integral);
    (*h)->SetFillColor(*c);
    mc_stack.Add(*h);
  }

  TCanvas *cv = new TCanvas();
  TPad *pad1 = new TPad("plotPad","",0.005,0.21,0.995,0.995);
  pad1->cd();
  //cv->Divide(1,2);
  //cv->cd(1);
  

  if(data_total->GetMaximum() > mc_stack.GetMaximum()) data_total->SetAxisRange(1e-2,data_total->GetMaximum()*1.2,"Y");
  else data_total->SetAxisRange(1e-2,mc_stack.GetMaximum()*1.2,"Y");
  data_total->SetMarkerStyle(8);
  data_total->Draw("PE1");
  mc_stack.Draw("HISTSAME");
  data_total->Draw("PE1SAME");
  
  TLegend leg(0.6,0.7,0.85,0.9);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(data_total,"Data","P");
  leg.AddEntry(mc[0],"Real Photons","F");
  leg.AddEntry(mc[1],"Real Electrons","F");
  leg.AddEntry(mc[2],"Fakes","F");
  leg.Draw("SAME");

  TLatex lbl(0.60,0.96,label);
  lbl.SetNDC();
  lbl.SetTextSize(0.045);
  lbl.SetTextColor(kBlack);

  TLatex var(0.12,0.96,xName);
  var.SetNDC();
  var.SetTextSize(0.045);
  var.SetTextColor(kBlack);

  lbl.Draw();
  var.Draw();

  
  TPad *pad2 = new TPad("ratioPad","",0.005,0.005,0.995,0.25);
  pad2->cd();
  //cv->cd(2);
  
  TH1F* ratio = (TH1F*)mc[0]->Clone("ratio");
  ratio->SetFillColor(0);
  
  for(int j=0;j<ratio->GetNbinsX();j++){
    float dataN = data_total->GetBinContent(j);
    float mcN = mc_total->GetBinContent(j);
    float mcE = mc_total->GetBinError(j);
    if(mcN){
      ratio->SetBinContent(j,dataN/mcN);
      if(dataN) ratio->SetBinError(j,dataN/mcN*sqrt(1/dataN+mcE*mcE/mcN/mcN));
      else ratio->SetBinError(j,0.6);
    }else{
      ratio->SetBinContent(j,1);
      ratio->SetBinError(j,0.6);
    }
  }
  
  ratio->SetYTitle("Data/MC");
  ratio->SetAxisRange(-0.1,2.0,"Y");
  ratio->SetXTitle(xName);
  ratio->SetFillColor(0);
  ratio->Draw("E1");
  
  cv->cd();
  pad1->Draw();
  pad2->Draw();


  return cv;  
}
