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



std::map<int,TString> getCatNames() {
  std::map<int,TString> catNames = { {0,"HighPt"}, {1,"Hbb"}, {2,"Zbb"}, {3,"HighRes"}, {4, "LowRes"}};
  return catNames;
}

int getBins(int reg, int* edges,bool isHT=false) {
  //edges = new int[4];



  //NOMINAL BINNING
  int HighPt_bins[] = {1,2,4,8,16};
  int Hbb_bins[] = {1,4};
  int Zbb_bins[] = {1,7};
  //int HighRes_bins[] = {1,4,6,15};
  int HighRes_bins[] = {1,3,6,15};
  int LowRes_bins[] = {1,2,3,6,14};


  /*
  //ALTERNATE BINNING
  int HighPt_bins[] = {1,2,4,8,16};
  int Hbb_bins[] = {1,5};
  int Zbb_bins[] = {1,7};
  //int HighRes_bins[] = {1,4,6,15};
  int HighRes_bins[] = {1,3,6,15};
  int LowRes_bins[] = {1,2,3,6,15};
  */
  /*
  //DEFAULT CORRECTION BINNING
  int HighPt_bins[] = {1,2,4,8,16};
  int Hbb_bins[] = {1,4};
  int Zbb_bins[] = {1,7};
  int HighRes_bins[] = {1,3,5,14};
  //int HighRes_bins[] = {1,2,4,13};
  int LowRes_bins[] = {1,2,3,6,14};
  */
  /*
  //V3 binning!
  int HighPt_bins[] = {1,3,8,16};
  int Hbb_bins[] = {1,4};
  int Zbb_bins[] = {1,4};
  int HighRes_bins[] = {1,2,5,14};
  int LowRes_bins[] = {1,3,6,14};
  */
  /*
    //Minus 50 GeV Binning
  int HighPt_bins[] = {1,2,7,15};
  int Hbb_bins[] = {1,3};
  int Zbb_bins[] = {1,3};
  int HighRes_bins[] = {1,2,4,13};
  int LowRes_bins[] = {1,2,5,13};

  */
  int HT_HighPt_bins[] = {1,20,25,30};
  int HT_Hbb_bins[] = {1,8};
  int HT_Zbb_bins[] = {1,11};
  int HT_HighRes_bins[] = {1,16,21,28};
  int HT_LowRes_bins[]  = {1,17,22};

  switch(reg) {
  case 0:
    std::copy(HighPt_bins,HighPt_bins+5,edges);
    if(isHT) std::copy(HT_HighPt_bins,HT_HighPt_bins+4,edges);
    //edges = HighPt_bins;
    return 5;
  case 1:
    std::copy(Hbb_bins,Hbb_bins+2,edges);
    if(isHT) std::copy(HT_Hbb_bins,HT_Hbb_bins+2,edges);
    //edges = Hbb_bins;
    return 2;
  case 2:
    std::copy(Zbb_bins,Zbb_bins+2,edges);
    if(isHT) std::copy(HT_Zbb_bins,HT_Zbb_bins+2,edges);
    //edges = Zbb_bins;
    return 2;
  case 3:
    std::copy(HighRes_bins,HighRes_bins+4,edges);
    if(isHT) std::copy(HT_HighRes_bins,HT_HighRes_bins+4,edges);
    //edges = HighRes_bins;
    return 4;
  case 4:
    std::copy(LowRes_bins,LowRes_bins+5,edges);
    if(isHT) std::copy(HT_LowRes_bins,HT_LowRes_bins+4,edges);
    //edges = LowRes_bins;
    return 5;
  default:
    edges = NULL;
    return 0;
  }
  return 0;
}


std::vector<TString> getDataShapeSys() {
  std::vector<TString> sys = { "bkgShape", "fit" };
  return sys;
}
std::vector<TString> getMCShapeSys() {
  std::vector<TString> sys = { "btag", "jec", "phoE", "sigE" };
  return sys;
}

std::vector<std::pair<float,float> > getMCScaleSys(TString type) {
  std::vector<std::pair<float,float> > sys;
  sys.push_back(std::make_pair(0.975,1.025));
  sys.push_back(std::make_pair(0.95,1.05));
  sys.push_back(std::make_pair(0.74,1.28));

  if(type=="ggH") {
      sys.push_back(std::make_pair(0.918,1.076));
      sys.push_back(std::make_pair(0.930,1.076));
  }
  if(type=="ttH") {
      sys.push_back(std::make_pair(0.918,1.076));
      sys.push_back(std::make_pair(0.930,1.076));
  }
  if(type=="vbfH") {
      sys.push_back(std::make_pair(0.992,1.003));
      sys.push_back(std::make_pair(0.972,1.026));
  }
  if(type=="wzH") {
      sys.push_back(std::make_pair(0.982,1.021));
      sys.push_back(std::make_pair(0.958,1.042));
  }

  return sys;
}



TH2F* makeSignalRegionPlot(TH2F* inputHist,TString name,const int nMRbins,int* binEdges,bool printRanges,std::vector<float>* data=0,bool silent=false) {
  TH2F* output = (TH2F*)inputHist->Clone(name);
  output->SetXTitle("M_{R} [GeV]");
  output->SetYTitle("R^{2}");
  output->SetZTitle("Events");
  for(int iXbin=0;iXbin<nMRbins;iXbin++) {
    for(int iYbin=0;iYbin<nMRbins-iXbin;iYbin++) {
      float sum = 0;
      for(int x=1;x<output->GetNbinsX()+1;x++) {
	for(int y=0;y<output->GetNbinsY()+1;y++) {
	  if(x<binEdges[iXbin]) continue;
	  if(iXbin<nMRbins-1 && x>=binEdges[iXbin+1]) continue;
	  if(y<iYbin+1) continue;
	  if(iYbin<nMRbins-iXbin-1 && y>iYbin+1) continue;
	  sum+=output->GetBinContent(x,y);
	}
      }
      if(!silent) {
	if(printRanges){
	  std::cout << std::setw(4) << output->GetXaxis()->GetBinLowEdge(binEdges[iXbin]) << " - ";
	  if(iXbin<nMRbins-1) std::cout << std::setw(4) << output->GetXaxis()->GetBinLowEdge(binEdges[iXbin+1]);
	  else std::cout << std::setw(4) << "3000";
	  
	  std::cout << "& " << std::setw(4) << output->GetYaxis()->GetBinLowEdge(iYbin+1) << " - ";
	  if(iYbin<nMRbins-iXbin-1) std::cout << std::setw(4) << output->GetYaxis()->GetBinLowEdge(iYbin+2);
	  else std::cout << "1.00";
	}
	std::cout << " & " << std::setw(12) << sum << std::setprecision(4) << std::endl;
      }
      if(data) data->push_back(sum);

      for(int x=1;x<output->GetNbinsX()+1;x++) {
	for(int y=0;y<output->GetNbinsY()+1;y++) {
	  if(x<binEdges[iXbin]) continue;
	  if(iXbin<nMRbins-1 && x>=binEdges[iXbin+1]) continue;
	  if(y<iYbin+1) continue;
	  if(iYbin<nMRbins-iXbin-1 && y>iYbin+1) continue;
	  output->SetBinContent(x,y,sum);
	}
      }
    }
  }
  return output;
}

void makeSignalRegionTable(std::vector<TH2F*>& inputHist,const int nMRbins,int* binEdges,std::vector<float>& sf) {
  const TH2F* output = inputHist.at(0);
  assert(output!=0);


  for(int iXbin=0;iXbin<nMRbins;iXbin++) {
    for(int iYbin=0;iYbin<nMRbins-iXbin;iYbin++) {
      for(int iHist=0;iHist<inputHist.size();iHist++) {
	assert(inputHist[iHist] != 0);
	float sum = 0;
	for(int x=1;x<output->GetNbinsX()+1;x++) {
	  for(int y=0;y<output->GetNbinsY()+1;y++) {
	    if(x<binEdges[iXbin]) continue;
	    if(iXbin<nMRbins-1 && x>=binEdges[iXbin+1]) continue;
	    if(y<iYbin+1) continue;
	    if(iYbin<nMRbins-iXbin-1 && y>iYbin+1) continue;
	    sum+=inputHist.at(iHist)->GetBinContent(x,y);
	  }
	}
	if(iHist==0) {
	  printf("% 5.0f -",output->GetXaxis()->GetBinLowEdge(binEdges[iXbin]));
	  if(iXbin<nMRbins-1) printf("% 5.0f",output->GetXaxis()->GetBinLowEdge(binEdges[iXbin+1]));
	  else printf( " 3000");
	  printf(" & %0.2f - ",output->GetYaxis()->GetBinLowEdge(iYbin+1));
	  if(iYbin<nMRbins-iXbin-1) printf("%0.2f",output->GetYaxis()->GetBinLowEdge(iYbin+2));
	  else printf("1.00");

	  /*
	  std::cout << std::setw(4) << output->GetXaxis()->GetBinLowEdge(binEdges[iXbin]) << " - ";
	  if(iXbin<nMRbins-1) std::cout << std::setw(4) << output->GetXaxis()->GetBinLowEdge(binEdges[iXbin+1]);
	  else std::cout << std::setw(4) << "3000";
	  
	  std::cout << "& " << std::setw(4) << output->GetYaxis()->GetBinLowEdge(iYbin+1) << " - ";
	  if(iYbin<nMRbins-iXbin-1) std::cout << std::setw(4) << output->GetYaxis()->GetBinLowEdge(iYbin+2);
	  else std::cout << "1.00";
	  */
	}
	if(sf.at(iHist)>=0) printf(" & $ % 7.2f \\pm % 7.3f $ ", sum, sqrt(sum*sf.at(iHist)));
	else printf(" & $ % 7.2f $ ", sum);
	  
	//std::cout << " & $"  << std::setw(6) << sum << " \\pm " << sqrt(sum*sf.at(iHist)) << std::setprecision(4) << " $ ";
      

      }
      std::cout << " \\\\ " << std::endl;
    }
  }
}


TH2F* getHist(TFile **f, int nFiles,TString tag) {
  TH2F* sum = (TH2F*)f[0]->Get(tag);
  for(int i=1;i<nFiles;i++) {
    sum->Add( (TH2F*)f[i]->Get(tag) );
  }
  return sum;
}


std::pair<float,float> getScaleFactor(int cat) {
  switch(cat) {
  case 0:
    return std::make_pair(0.160,0.005);
  case 1:
    return std::make_pair(0.156,0.045);
  case 2:
    return std::make_pair(0.185,0.053);
  case 3:
    return std::make_pair(0.165,0.003);
  case 4:
    return std::make_pair(0.291,0.0034);


  }
}


void makeSigRegionTable(TFile **f, int nFiles,bool *blind,float *sfs,TString *infixes) {
  int binEdges[5];
  std::map<int,TString> catNames = getCatNames();
  for(int i=0;i<5;i++) {
    int nBins = getBins(i,binEdges);
    std::vector<TH2F*> hists;
    std::vector<float> scales;
    std::cout << std::endl << catNames[i] << std::endl;     
    for(int iFile=0;iFile<nFiles;iFile++) {
      TString suffix= (blind[iFile] ? "_SidebandRegion" : "_SignalRegion");
      //std::cout << "data_"+infixes[i]+catNames[i]+suffix << std::endl;
      hists.push_back(getHist(&f[iFile],1,"data_"+infixes[iFile]+catNames[i]+suffix));      
      scales.push_back( (sfs[iFile]==0 ? getScaleFactor(i).first : sfs[iFile]) );
    }   
    std::cout << hists.at(0) << std::endl;
    assert(hists.size()!=0);
    makeSignalRegionTable(hists,nBins,binEdges,scales);
  }
}

void makeYieldTable(TString dir,bool noSMS=false) {
  const int nFiles = 6;

  TFile *f[nFiles] = { TFile::Open(dir+"/data.root"),TFile::Open(dir+"/SMHiggs_SUM.root"),TFile::Open(dir+"/sms_ChiWH.root"),TFile::Open(dir+"/sms_ChiWH.root"),TFile::Open(dir+"/sms_ChiHH.root"),TFile::Open(dir+"/sms_ChiHH.root") };
  bool blind[nFiles] = {true,false,false,false,false,false};
  float  sfs[nFiles] = {0,992./(96290+99885+100320+93312),4.76*2.28E-3*19780/26789.,1.32*2.28E-3*19780/11612.,4.12*2.28E-3*19780/49898.,0.785*2.28E-3*19780/11618.};
  TString infixes[nFiles] = {"","","0_125_","0_200_","0_125_","0_200_"};

  makeSigRegionTable(f,(noSMS ? 2 : nFiles),blind,sfs,infixes);
}

void makeSMHiggsTable(TString dir) {
  const int nFiles = 5;

  TFile *f[nFiles] = { TFile::Open(dir+"/SMHiggs_SUM.root"),TFile::Open(dir+"/ggH.root"), TFile::Open(dir+"/vbfH.root"), TFile::Open(dir+"/wzH.root"), TFile::Open(dir+"/ttH.root")};
  bool blind[nFiles] = {false,false,false,false,false};
  float  sfs[nFiles] = {992./(96290+99885+100320+93312),2.28E-3*19780*19.78/96290,2.28E-3*19780*1.578/99885,2.28E-3*19780*0.7046/100320,2.28E-3*19780*0.1293/93312 };
  TString infixes[nFiles] = {"","","","",""};

  makeSigRegionTable(f,nFiles,blind,sfs,infixes);
}

void makeDataSMHiggsTable(TString dir) {
  const int nFiles = 6;

  TFile *f[nFiles] = { TFile::Open(dir+"/data.root"),TFile::Open(dir+"/data.root"),TFile::Open(dir+"/ggH.root"), TFile::Open(dir+"/vbfH.root"), TFile::Open(dir+"/wzH.root"), TFile::Open(dir+"/ttH.root")};
  bool blind[nFiles] = {false,true,false,false,false,false};
  float  sfs[nFiles] = {-1,0,2.28E-3*19780*19.78/96290,2.28E-3*19780*1.578/99885,2.28E-3*19780*0.7046/100320,2.28E-3*19780*0.1293/93312 };
  TString infixes[nFiles] = {"","","","","",""};

  makeSigRegionTable(f,nFiles,blind,sfs,infixes);
}


int makeSigRegionFromFile(TFile **f, int nFiles,TString outputDir,TString tag,TString histInfix="",float norm=1,bool blind=true,bool save=true,bool printRanges=true,bool make1D=false,bool isHT=false) {

  std::map<int,TString> catNames = getCatNames();
  int binEdges[10];
  TString suffix = (blind ? "_SidebandRegion" : "_SignalRegion");
  std::vector<float> *data = 0;
  TFile *ff = 0;
  if(save) ff = new TFile(outputDir+"/signalRegions_"+tag+suffix+".root","RECREATE");

  for(int iCat=0;iCat<5;iCat++) {
    int nBins = getBins(iCat,binEdges,isHT);
    if(make1D) { 
      delete data;
      data = new std::vector<float>;
    }
    std::cout << std::endl << catNames[iCat] <<std::endl;
    TH2F* hist = getHist(f,nFiles,"data_"+histInfix+catNames[iCat]+suffix);
    TH2F* hist_sig = makeSignalRegionPlot(hist,"sigregions_"+catNames[iCat]+suffix,nBins,binEdges,printRanges,data);
  
    hist_sig->Scale(norm);
    
    if(!save) continue;

    TCanvas cv;

    cv.SetLogz();

    hist_sig->SetAxisRange(0.001,1,"Y");
    hist_sig->Draw("COLZ");
    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catNames[iCat]+".png");
    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catNames[iCat]+".pdf");
    cv.SetLogy();
    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catNames[iCat]+"_LOG.png");
    cv.SaveAs(outputDir+"/signalRegions_"+tag+suffix+"_"+catNames[iCat]+"_LOG.pdf");


    ff->cd();
    hist_sig->Write();
  }

  ff->Close();

  return 0;
}

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
  */

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
	*/
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

