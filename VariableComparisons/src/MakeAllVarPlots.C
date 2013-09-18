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
#include "TMath.h"

#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <array>
#include <assert.h>
#include <string>
#include <tuple>

#include "../include/ArgParser.hh"
#include "include/HistogramStack.hh"


#include <cstdlib>

//local includes
#include "include/weightManager.hh"
#include "include/varCorrector.hh"
#include "include/plotManager.hh"
#include "src/getCategories.C"


void setstyle();

TCanvas *makeCanvas(std::array<TH1F*,3> data,std::array<TH1F*,3> mc, TString xName,TString label="");

enum RunPeriod {kAll,kAB,kABC,kD};
struct options_t {
  bool isTest = false;
  bool useSherpa = false;
  RunPeriod period = kAll;
};

void writeRootFile(TString fileName,options_t opt);
void DrawFromRootFile(TString fileName,options_t opt); //bool useSherpa,bool test);



int main(int argc,char** argv){
  ArgParser a(argc,argv);


  a.addArgument("mode",ArgParser::required,"specify the operation mode (save,draw)");
  a.addArgument("rootFileName",ArgParser::required,"specify the name of the root file");
  a.addLongOption("sherpa",ArgParser::noArg,"use SHERPA diphoton jets sample");
  a.addLongOption("test",ArgParser::noArg,"Only process 1 MC file (useful for testing)");
  a.addLongOption("run",ArgParser::reqArg,"Specify the run period to plot [ABC,D]");
  a.addShortOption('h',ArgParser::noArg,"print help");
  std::string ret;
  if(a.process(ret)!=0) {
    std::cout << "\nInvalid Option " << ret << std::endl << std::endl;;
    a.printOptions(argv[0]);
    return 0;
  }

  if(a.shortFlagPres('h')) {
    a.printOptions(argv[0]);
    return 0;
  }
  
  enum OpMode{kNone,kSave,kDraw};
  OpMode mode = kNone;

  if(a.getArgument("mode").compare("save")==0) mode=kSave;
  if(a.getArgument("mode").compare("draw")==0) mode=kDraw;
  
  if(mode==kNone) {
    std::cout << "\nInvalid Mode Speficied:  " << a.getArgument("mode") << std::endl << std::endl;
    return -1;
  }
  TString fileName = a.getArgument("rootFileName");

  options_t opt;
  opt.useSherpa = a.longFlagPres("sherpa");
  opt.isTest    = a.longFlagPres("test");
  opt.period = kAll;
  if(a.longFlagPres("run")) {
    if(a.getLongFlag("run").compare("AB")==0) opt.period = kAB;
    if(a.getLongFlag("run").compare("ABC")==0) opt.period = kABC;
    if(a.getLongFlag("run").compare("D")==0) opt.period = kD;    
  }

  switch(mode) {
  case kSave:
    writeRootFile(fileName,opt); //useSherpa,isTest);
    break;
  case kDraw:
    DrawFromRootFile(fileName,opt); //useSherpa,isTest);
    break;
  }


  return 0;
}

void writeRootFile(TString fileName, options_t opt) { //bool useSherpa,bool test) {
  bool useSherpa = opt.useSherpa;
  bool test    = opt.isTest;

  if(useSherpa) std::cout << "SHERPA!!" << std::endl;


  float lumi=6.3;
  weightManager weights;

  plotManager mcPlotter("mc");
  plotManager dataPlotter("data");
  mcPlotter.setUse4Cat();
  dataPlotter.setUse4Cat();

  setstyle();

  std::vector<TString> samples = {
    "QCD_Pt30to40.root",
    "QCD_Pt40.root",
    "GJets_Pt20to40.root",
    "GJets_Pt40.root",
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

  int Nsamples = samples.size();

  std::vector<TString> data;
  if(opt.period == kABC || opt.period == kAB || opt.period == kAll) {
    data.push_back("Photon_Run2012A_22Jan2013.root");
    data.push_back("DoublePhoton_Run2012B_22Jan2013.root");
  }
  if(opt.period ==kABC || opt.period ==kAll) {
    data.push_back("DoublePhoton_Run2012C_22Jan2013.root");
  }
  if(opt.period == kD || opt.period == kAll) {
    data.push_back("DoublePhoton_Run2012D_22Jan2013.root");
  }


  int Ndata = data.size();


  std::vector<TString> globalVetos;
  globalVetos.push_back("etaSC > -1.775 && etaSC < -1.764 && phi > 1.2 && phi < 1.6");
  globalVetos.push_back("etaSC > 1.591 && etaSC < 1.595 && phi > -2.06 && phi < -2.045");
  //globalVetos.push_back("etaSC > -0.9 && etaSC < -0.7 && phi > 2.9 && phi < 3.1");
  globalVetos.push_back("etaSC > 1.74 && etaSC < 1.76 && phi > 2.1 && phi < 2.15");
  globalVetos.push_back("etaSC > 1.564 && etaSC < 1.565 && phi > 0.528 && phi < 0.532");

  //output->Draw("phi:etaSC","!(etaSC > -1.8 && etaSC < -1.76 && phi > 1.2 && phi < 1.5) && !(etaSC>1.591 && etaSC<1.595 && phi > -2.06 && phi < -2.045) && !(etaSC > 1.74 && etaSC < 1.76 && phi > 2.1 && phi < 2.15) && !(etaSC > 1.564 && etaSC < 1.565 && phi > 0.528 && phi < 0.532) && !(etaSC > -1.6888 && etaSC < -1.6868 && phi > 1.8 && phi < 1.802) && sieie < 1e-5 && abs(etaSC) > 1.557 && abs(etaSC)<2

  auto catInfo = getCategories();
  auto cats = catInfo.first;
  auto catNames = catInfo.second;



  assert(cats.size() == catNames.size());

  int Ncat = cats.size();

  if(test) {
    Ncat = Ndata = Nsamples = 1;
  }

  for(int i=0;i<Ncat;i++){
    mcPlotter.addCategory(catNames[i],cats[i]);
    dataPlotter.addCategory(catNames[i],cats[i]);
  }

  mcPlotter.addVetos(&globalVetos);
  dataPlotter.addVetos(&globalVetos);

  auto vars = getVars();
  for(auto it : vars) {
    mcPlotter.addVariable  (translateVar(std::get<0>(it)),std::get<0>(it),std::get<1>(it),std::get<2>(it),std::get<3>(it));
    dataPlotter.addVariable  (translateVar(std::get<0>(it)),std::get<0>(it),std::get<1>(it),std::get<2>(it),std::get<3>(it));
  }

  int Nvar = vars.size();

  
  TFile *data_pu_file = new TFile("/home/amott/HggApp/PU/truepileup_69400mb_2012_RunABCD_22Jan2013.root");
  TH1D  *data_pu_hist = (TH1D*)data_pu_file->Get("pileup_hist")->Clone("data_pu_hist");
  data_pu_hist->SetDirectory(0);
  data_pu_file->Close();
  std::cout << data_pu_hist->GetName() <<std::endl;

  mcPlotter.setTargetPU(data_pu_hist);
  

  TFile outputFile(fileName,"RECREATE");
  outputFile.cd();
  mcPlotter.buildHistograms();
  for(int iSample=0;iSample<Nsamples;iSample++){
    TString sampleName = samples[iSample];
    sampleName.Remove(sampleName.Last('.'));
    
    TFile mc_pu_file("/home/amott/HggApp/PU/"+samples[iSample]);
    TH1D* mc_pu_hist = (TH1D*)mc_pu_file.Get("pileup_hist")->Clone("mc_pu_hist");
    mc_pu_hist->SetDirectory(0);
    mc_pu_file.Close();
    mcPlotter.setMCPU(mc_pu_hist);
    
    TFile *f = new TFile("output/"+samples[iSample]);
    TChain *fChain = (TChain*)f->Get("output");
    mcPlotter.processChain(fChain,weights.getWeight(sampleName,lumi));
    f->Close();
    mcPlotter.setMCPU(0);
    delete mc_pu_hist;
  }
  delete data_pu_hist;
  mcPlotter.setTargetPU(0);

  outputFile.cd();
  dataPlotter.buildHistograms();
  for(int iData=0;iData<Ndata;iData++){
    TFile *f = new TFile("output/"+data[iData]);
    TChain *fChain = (TChain*)f->Get("output");
    dataPlotter.processChain(fChain,1);
  }


//   for(int iCat=0;iCat<Ncat;iCat++){
//     for(int iVar=0;iVar<Nvar;iVar++){
//       std::array<TH1F*,3> mc = mcPlotter.getHistogram(catNames[iCat],std::get<0>(vars[iVar]));
//       //std::array<TH1F*,3> mcCorr = mcPlotter.getHistogram(catNames[iCat],"corr_"+vars[iVar]);
//       std::array<TH1F*,3> data = dataPlotter.getHistogram(catNames[iCat],std::get<0>(vars[iVar]));
//       //std::array<TH1F*,3> dataCorr = dataPlotter.getHistogram(catNames[iCat],"corr_"+vars[iVar]);
//       for(auto h : mc)       if(h) {h->Write();}
//       //for(auto h : mcCorr)   h->Write();
//       for(auto h : data)     if(h) {h->Write();}
//       //for(auto h : dataCorr) h->Write();
//     }
//   }
  mcPlotter.saveAll(&outputFile);
  dataPlotter.saveAll(&outputFile);
  outputFile.Close();
} //void writeRootFile


void DrawFromRootFile(TString fileName,options_t opt) { //bool useSherpa,bool isTest) {
  auto catInfo = getCategories();
  auto vars     = getVars();

  bool useSherpa = opt.useSherpa;
  bool isTest    = opt.isTest;

  setstyle();

  std::vector<TString> catNames = catInfo.second;
  int Ncat = catNames.size();
  int Nvar = vars.size();

  TFile *file = new TFile(fileName);

  for(auto catNameIt : catNames) {
    for(auto varIt : vars) {
      TString varName = std::get<0>(varIt);
      std::array<TH1F*,3> mc,data;


    }//for(auto varIt : vars)
  }//  for(auto catNameIt : catNames)

  const int nTypes=3;

  TString folder  = "figs";
  if(opt.period == kAB) folder = folder+"_RunAB";
  if(opt.period == kABC) folder = folder+"_RunABC";
  if(opt.period == kD)   folder = folder+"_RunD";
  folder = folder+"/";
  if(useSherpa) folder  = folder+"SHERPA/";


  std::array<TString,nTypes> histSuffix = {"realPho","realEle","fake"};
  for(auto catIt : catNames) {
    for(auto varIt : vars) {
      std::array<TH1F*,nTypes> mc,data;
      TString transVar = translateVar(std::get<0>(varIt));

      for(int i=0;i<nTypes;i++) {
	mc[i] = (TH1F*)file->Get( transVar+"_"+catIt+"_mc_"+histSuffix[i] );
	data[i] = (TH1F*)file->Get( transVar+"_"+catIt+"_data_"+histSuffix[i] );
      }
      if(mc[0]==0) continue;
      TCanvas *cv = makeCanvas(data,mc,std::get<0>(varIt),catIt);


      TString saveVar = transVar;
      if(isTest) saveVar+="__TEST";
      
      if(useSherpa) {
	saveVar = "SHERPA__"+saveVar;
      }

      TString path = folder+saveVar+"_"+catIt;

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
    }
  }  
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
  std::cout << data[0] << std::endl;
  TH1F *data_total = (TH1F*)data[0]->Clone("data_total");
  data_total->Add(data[1]);
  data_total->Add(data[2]);
  
  float data_norm = data_total->Integral();
  std::cout << "# data events:  " << data_norm << std::endl;
  
  
  HistogramStack<TH1F> mc_stack;
  //compute the total integral for normalization
  double mc_integral = 0;
  for(auto h : mc) mc_integral+=h->Integral();

  //make the stack
  std::array<Color_t,3> colors = {kBlue,kGreen,kRed};
  assert(colors.size() == mc.size());

  auto h = mc.begin();
  auto c = colors.begin();

  for(; h!=mc.end(); h++,c++) {
    //scale MC to data
    (*h)->Scale(data_norm/mc_integral);
    (*h)->SetFillColor(*c);
    mc_stack.Add(**h);
  }

  TCanvas *cv = new TCanvas();
  TPad *pad1 = new TPad("plotPad","",0.005,0.21,0.995,0.995);
  pad1->cd();
  //cv->Divide(1,2);
  //cv->cd(1);
  
  double KS = mc_stack.getTotal()->KolmogorovTest(data_total,"");
  KS = TMath::Log(KS);
  std::cout << "KS: " << KS << std::endl;


  //if(data_total->GetMaximum() > mc_stack.GetMaximum()) data_total->SetAxisRange(1e-2,data_total->GetMaximum()*1.2,"Y");
  //else data_total->SetAxisRange(1e-2,mc_stack.GetMaximum()*1.2,"Y");
  data_total->SetMarkerStyle(8);
  data_total->Draw("PE1");
  mc_stack.Draw("HISTSAME");
  data_total->Draw("PE1SAME");

  //cv->cd();
  TLatex *lbl = new TLatex(0.60,0.96,label);
  lbl->SetNDC();
  lbl->SetTextSize(0.045);
  lbl->SetTextColor(kBlack);

  TLatex *var = new TLatex(0.12,0.96,xName);
  var->SetNDC();
  var->SetTextSize(0.045);
  var->SetTextColor(kBlack);

  TLatex *KSval = new TLatex(0.12,0.91, Form("Log(KS Prob): %0.3f",KS));
  KSval->SetNDC();
  KSval->SetTextSize(0.045);
  KSval->SetTextColor(kBlack);

  lbl->Draw();
  var->Draw();
  KSval->Draw();

  TLegend *leg = new TLegend(0.6,0.7,0.85,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(data_total,"Data","P");
  leg->AddEntry(mc[0],"Real Photons","F");
  leg->AddEntry(mc[1],"Real Electrons","F");
  leg->AddEntry(mc[2],"Fakes","F");
  leg->Draw("SAME");


  
  TPad *pad2 = new TPad("ratioPad","",0.005,0.005,0.995,0.25);
  pad2->cd();
  //cv->cd(2);
  
  TH1F* ratio = (TH1F*)mc[0]->Clone("ratio");
  ratio->SetFillColor(0);
  
  for(int j=0;j<ratio->GetNbinsX();j++){
    float dataN = data_total->GetBinContent(j);
    float mcN = mc_stack.getTotal()->GetBinContent(j);
    float mcE = mc_stack.getTotal()->GetBinError(j);
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
