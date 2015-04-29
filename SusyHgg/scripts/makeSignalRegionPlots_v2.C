#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TColor.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include "TRandom3.h"

#include <map>
#include <cmath>
#include "assert.h"

#include "SigRegionBinning.h"





int makeSigRegionFromFile(TFile **f, int nFiles,TString outputDir,TString tag,TString histInfix="",bool blind=true) {

  TString suffix = (blind ? "_SidebandRegion" : "_SignalRegion");
  for(int iCat=0;iCat<5;iCat++) {
    TString catName = SigRegionBinning::getRegionName(SigRegionBinning::br(iCat));
    TH2F* hist_sig = SigRegionBinning::makeHistogram(SigRegionBinning::br(iCat),"sigregions_"+catName+suffix);

    for(int iFile=0; iFile<nFiles; iFile++) {
      SigRegionBinning::addToSignalRegionHistogram(hist_sig, (TH2F*) f[iFile]->Get("data_"+histInfix+catName+suffix));
    }
    
    SigRegionBinning::formatSigRegionPlot(hist_sig);
    TCanvas cv;
    cv.SetLogz();

    hist_sig->SetAxisRange(0.001,1,"Y");
    hist_sig->Draw("COLZ");
    SigRegionBinning::drawSigRegionLines(&cv,SigRegionBinning::br(iCat));

    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catName+".png");
    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catName+".pdf");
    cv.SetLogy();
    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catName+"_LOG.png");
    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catName+"_LOG.pdf");
  }

  return 0;
}

int makeSigRegionFromFile(TFile *f, TString outputDir,TString tag,TString histInfix="",bool blind=true) {
  TFile *ff[] ={f};
  return makeSigRegionFromFile(ff,1,outputDir,tag,histInfix,blind);
}

int makeTotalBkg(TString dir="./",bool percent=false) {
  TFile dataFile(dir+"/data.root");
  TFile smhFile(dir+"/SMHiggs_SUM.root");

  if(!smhFile.IsOpen()) { 
    std::cout << "Cannot open SMHiggs_SUM.root (does it exist?)" << std::endl;
    return -1;
  }

  Double_t Red[] = {0.00, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00};
  Double_t Green[] ={0.00, 0.70, 0.90, 1.00, 0.90, 0.70, 0.00};
  Double_t Blue[] = {1.00, 1.00, 1.00, 1.00, 0.90, 0.70, 0.00};
  Double_t Length[] =  {0.00, 0.20, 0.35, 0.50, 0.65, 0.8, 1.00};


  for(int iCat=0;iCat<5;iCat++) {
    TString catName = SigRegionBinning::getRegionName(SigRegionBinning::br(iCat));
    TH2F* obs_hist_sig = SigRegionBinning::makeHistogram(SigRegionBinning::br(iCat),"obs_"+catName, (TH2F*) dataFile.Get("data_"+catName+"_SignalRegion"));
    TH2F* bkg_hist_sig = SigRegionBinning::makeHistogram(SigRegionBinning::br(iCat),"bkg_"+catName, (TH2F*) dataFile.Get("data_"+catName+"_SidebandRegion"));
    TH2F* smh_hist_sig = SigRegionBinning::makeHistogram(SigRegionBinning::br(iCat),"smh_"+catName, (TH2F*) smhFile.Get("data_"+catName+"_SignalRegion"));
        
    SigRegionBinning::formatSigRegionPlot(obs_hist_sig);
    SigRegionBinning::formatSigRegionPlot(bkg_hist_sig);
    SigRegionBinning::formatSigRegionPlot(smh_hist_sig);

    TH2F* nsig = (TH2F*)bkg_hist_sig->Clone("nsig");
    nsig->Add(smh_hist_sig);

    TH2F* perc = (TH2F*)bkg_hist_sig->Clone("perc");
    perc->Add(smh_hist_sig);

    for(int iXbin=1; iXbin<=nsig->GetNbinsX(); iXbin++) {
      for(int iYbin=1; iYbin<=nsig->GetNbinsY(); iYbin++) {
	float obs = obs_hist_sig->GetBinContent(iXbin,iYbin);
	float bkg = nsig->GetBinContent(iXbin,iYbin);
	nsig->SetBinContent(iXbin,iYbin, (obs-bkg)/sqrt(bkg));
			    //TMath::NormQuantile(TMath::Gamma(obs+1,bkg)/TMath::Factorial(obs)));
	perc->SetBinContent(iXbin,iYbin, int(100*(obs-bkg)/bkg)/100.);	
      }
    }
    TCanvas cv;
    //cv.SetLogz();

    if(percent) {
      perc->Draw("COLZ TEXT45");
      SigRegionBinning::drawSigRegionLines(&cv,SigRegionBinning::br(iCat));

      cv.SaveAs(dir+"/figs/signalRegions_"+catName+"_PERDIFF.png");
      cv.SaveAs(dir+"/figs/signalRegions_"+catName+"_PERDIFF.pdf");
    } else {

      TColor::CreateGradientColorTable(7,Length,Red,Green,Blue,999);
      nsig->SetMinimum(-5.1);
      nsig->SetMaximum(5.1);
      nsig->SetContour(999);
      //nsig->SetAxisRange(0.001,1,"Y");
      nsig->Draw("COLZ");
      SigRegionBinning::drawSigRegionLines(&cv,SigRegionBinning::br(iCat));
      
      cv.SaveAs(dir+"/figs/signalRegions_"+catName+"_NSIG.png");
      cv.SaveAs(dir+"/figs/signalRegions_"+catName+"_NSIG.pdf");

    }

  }
}

/*
std::pair<TH2F*,TH2F*> sigRegionCompBox(TFile* f_data_obs,TFile* f_data_bkg,TFile* f_smtot,TString box){
  TH2F* data_obs = (TH2F*)f_data_obs->Get( Form("sigregions_%s_SignalRegion",box.Data()) );
  TH2F* data_bkg = (TH2F*)f_data_bkg->Get( Form("sigregions_%s_SidebandRegion",box.Data()) );
  TH2F* smtot = (TH2F*)f_smtot->Get( Form("sigregions_%s_SignalRegion",box.Data()) );

  TH2F* comp = (TH2F*)data_obs->Clone("data_diff_"+box);
  TH2F* nsig = (TH2F*)data_obs->Clone("data_nsig_"+box);

  comp->SetXTitle("M_{R} [GeV]");
  nsig->SetXTitle("M_{R} [GeV]");
  comp->SetYTitle("R^{2}");
  nsig->SetYTitle("R^{2}");
  comp->SetZTitle("Events");
  nsig->SetZTitle("Events");

  
  for(int iX=1;iX<comp->GetNbinsX()+1;iX++) {
    for(int iY=1;iY<comp->GetNbinsY()+1;iY++) {
      float v_data_obs = data_obs->GetBinContent(iX,iY);
      float v_data_bkg = data_bkg->GetBinContent(iX,iY);
      float v_smtot = smtot->GetBinContent(iX,iY);
      comp->SetBinContent(iX,iY, (v_data_obs-v_data_bkg-v_smtot));
      nsig->SetBinContent(iX,iY, (v_data_obs-v_data_bkg-v_smtot)/sqrt(v_data_bkg+v_smtot));
    }
  }
    
  return std::make_pair(comp,nsig);
}

void makeSigRegionComparison(TFile* data_obs,TFile* data_bkg,TFile* smtot,TString outputDir) {
  /*
 KEY: TH2F     sigregions_HighPt_SignalRegion;1
 KEY: TH2F     sigregions_Hbb_SignalRegion;1
 KEY: TH2F     sigregions_Zbb_SignalRegion;1
 KEY: TH2F     sigregions_HighRes_SignalRegion;1
 KEY: TH2F     sigregions_LowRes_SignalRegion;1


  Double_t Red[] = {0.00, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00};
  Double_t Green[] ={0.00, 0.70, 0.90, 1.00, 0.90, 0.70, 0.00};
  Double_t Blue[] = {1.00, 1.00, 1.00, 1.00, 0.90, 0.70, 0.00};
  Double_t Length[] =  {0.00, 0.20, 0.35, 0.50, 0.65, 0.8, 1.00};


  const int nBox=5;
  //const int nBox=4;
  TString boxes[nBox] = { "HighPt","Hbb","Zbb","HighRes","LowRes" };

  TCanvas cv;
for(int i=0;i<nBox;i++) {
    std::pair<TH2F*,TH2F*> hist = sigRegionCompBox(data_obs,data_bkg,smtot,boxes[i]);
    hist.first->SetAxisRange(0.001,1,"Y");
    hist.first->Draw("COLZ");
    cv.SetLogy(0);
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp.png");
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp.pdf");
    cv.SetLogy(1);
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp_LOG.png");
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp_LOG.pdf");
  }

  TColor::CreateGradientColorTable(7,Length,Red,Green,Blue,999);
  for(int i=0;i<nBox;i++) {
    std::pair<TH2F*,TH2F*> hist = sigRegionCompBox(data_obs,data_bkg,smtot,boxes[i]);
    hist.second->SetMinimum(-5.1);
    hist.second->SetMaximum(5.1);
    hist.second->SetContour(999);
    hist.second->SetAxisRange(0.001,1,"Y");
    hist.second->Draw("COLZ");
    cv.SetLogy(0);
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp_NSIG.png");
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp_NSIG.pdf");
    cv.SetLogy(1);
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp_NSIG_LOG.png");
    cv.SaveAs(outputDir+"/signalRegions_"+boxes[i]+"_obs_minus_exp_NSIG_LOG.pdf");

  }
}


float getPVal(float exp, float sf, float sfE,float obs) {
  TRandom3 rng(0);

  size_t count=0;

  float diff = fabs(obs-exp);

  size_t nToy = 100;
  while(count<100) {
    nToy*=10;
    count=0;
    for(size_t i=0;i<nToy;i++) {
      float thisExp = rng.Poisson(exp/sf);
      float thisSF  = rng.Gaus(sf,sfE);
      thisExp*=thisSF;
      float r = rng.Poisson(thisExp);

      //if( rng.Poisson(thisExp) >=obs ) count++;
      if( fabs(r-exp) >= diff ) count++;
    }
  }
  return float(count)/nToy;
}

std::vector<float> getSystematicsErrors(std::vector<float>* nominal,std::vector<std::vector<float>*>* errors) {
  std::vector<float> systs(nominal->size(),0);

  for(int iSyst=0;iSyst<errors->size();iSyst++) {
    if(errors->at(iSyst)->size()!=nominal->size()) {
      assert(false);
    }
    for(int iBin=0;iBin<nominal->size();iBin++) {
      systs.at(iBin)+=TMath::Power(nominal->at(iBin)-errors->at(iSyst)->at(iBin),2);
    }
  }
  return systs;
}

std::vector<std::vector<float> >* getSysts(TFile** files,int nFiles,TString name,bool isMC) {
  std::vector<TString> shapeSysts;
  if(isMC) shapeSysts = getMCShapeSys();
  else shapeSysts = getDataShapeSys();

  
}

void makeAll(TString folder,bool makePlots=true,bool isHT=false,float scale_ttH=1) {
  TFile *dataFile = new TFile(folder+"/data.root");
  TFile *data[] = {dataFile};

  if(makePlots) makeSigRegionFromFile(data,1,folder+"/figs/","DATA","",1,true,true,false,false,isHT);
  if(makePlots)   makeSigRegionFromFile(data,1,folder+"/figs/","DATA","",1,false,true,false,false,isHT);

  TFile *ggHFile  = new TFile(folder+"/ggH.root");
  TFile *vbfHFile = new TFile(folder+"/vbfH.root");
  TFile *wzHFile  = new TFile(folder+"/wzH.root");
  TFile *ttHFile  = new TFile(folder+"/ttH.root");
  TFile *sm [] = {ggHFile,vbfHFile,wzHFile,ttHFile};
  if(makePlots)   makeSigRegionFromFile(sm,4,folder+"/figs/","SMTot","",1,false,true,false,false,isHT);



  int binEdges[5];
  std::map<int,TString> catNames = getCatNames();
  for(int i=0;i<5;i++) {
    int nBins = getBins(i,binEdges,isHT);    
    TH2F* h_data_obs = getHist(data,1,"data_"+catNames[i]+"_SignalRegion");      
    std::vector<float> v_data_obs;
    TH2F* tmp = makeSignalRegionPlot(h_data_obs,"tmp",nBins,binEdges,false,&v_data_obs,true);
    delete tmp;

    TH2F* h_data_bkg = getHist(data,1,"data_"+catNames[i]+"_SidebandRegion");      
    std::vector<float> v_data_bkg;
    tmp = makeSignalRegionPlot(h_data_bkg,"tmp",nBins,binEdges,false,&v_data_bkg,true);
    delete tmp;

    TH2F* h_smtot = getHist(sm,3,"data_"+catNames[i]+"_SignalRegion");     
    TH2F* h_ttH = getHist(&sm[3],1,"data_"+catNames[i]+"_SignalRegion");     
    h_ttH->Scale(scale_ttH);
    h_smtot->Add(h_ttH);

    std::vector<float> v_smtot;
    tmp = makeSignalRegionPlot(h_smtot,"tmp",nBins,binEdges,false,&v_smtot,true);
    

    int j=-1;
    for(int iXbin=0;iXbin<nBins;iXbin++) {
      for(int iYbin=0;iYbin<nBins-iXbin;iYbin++) {
	++j;
	printf("% 4.0f -",tmp->GetXaxis()->GetBinLowEdge(binEdges[iXbin]));
	if(iXbin<nBins-1) printf("% 4.0f",tmp->GetXaxis()->GetBinLowEdge(binEdges[iXbin+1]));
	else printf("3000");
	//std::cout << std::fixed << tmp->GetXaxis()->GetBinLowEdge(binEdges[iXbin]) << " - ";
	//if(iXbin<nBins-1) std::cout << std::setw(5) << std::fixed << tmp->GetXaxis()->GetBinLowEdge(binEdges[iXbin+1]);
	//else std::cout << std::setw(5) << std::fixed << "3000";
	  
	printf(" & %0.2f - ",tmp->GetYaxis()->GetBinLowEdge(iYbin+1));

	if(iYbin<nBins-iXbin-1) printf("%0.2f",tmp->GetYaxis()->GetBinLowEdge(iYbin+2));
	else printf("1.00");
	//std::cout << "& " << std::setw(5) << std::fixed << tmp->GetYaxis()->GetBinLowEdge(iYbin+1) << " - ";
	//if(iYbin<nBins-iXbin-1) std::cout << std::setw(5) << std::fixed << tmp->GetYaxis()->GetBinLowEdge(iYbin+2);
	//else std::cout << "1.00";

	float obs = v_data_obs.at(j);

	float sideband = v_data_bkg.at(j);
	float sidebandE = sqrt(v_data_bkg.at(j) * getScaleFactor(i).first);
	float sm = v_smtot.at(j);
	float smE = sqrt(v_smtot.at(j) * 992./(96290+99885+100320+93312));

	float exp = sideband+sm;
	float expE = sqrt(sidebandE*sidebandE+smE*smE);

	float pv = getPVal(exp,getScaleFactor(i).first,getScaleFactor(i).second,obs);
	float nsig = 0;
	if(pv>0 && pv <1) nsig = TMath::NormQuantile(pv/2.);
	if(obs > exp) nsig*=-1;

	printf(" & % 4.0f & $ % 4.2f \\pm %4.2f $ & %0.5f & %0.2f \\\\ \n",obs, exp,expE,pv,nsig);
	/*
	std::cout << " & " << std::setw(12) << obs << std::setprecision(1);
	std::cout << " & $" << std::setw(12) << exp << " \\pm " << expE << " $ " << std::setprecision(3);
	std::cout << " & " << pv << std::setprecision(3);
	std::cout << " & " << nsig << std::setprecision(3);
	std::cout << " \\\\" << std::endl;

      }
    }
    std::cout << std::endl;
    delete tmp;
  }

  dataFile->Close();
  ggHFile->Close();
  vbfHFile->Close();
  wzHFile->Close();
  ttHFile->Close();


  TFile *data_obs = new TFile(folder+"/figs/signalRegions_DATA_SignalRegion.root");
  TFile *data_bkg = new TFile(folder+"/figs/signalRegions_DATA_SidebandRegion.root");
  TFile *smtot    = new TFile(folder+"/figs/signalRegions_SMTot_SignalRegion.root");


  //void makeSigRegionComparison(TFile* data_obs,TFile* data_bkg,TFile* smtot,TString outputDir) 
  if(makePlots)  makeSigRegionComparison(data_obs,data_bkg,smtot,folder+"/figs/");

  data_obs->Close();
  data_bkg->Close();
  smtot->Close();
}

*/
