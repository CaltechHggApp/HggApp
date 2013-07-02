#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooGlobalFunc.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include "RooStats/SPlot.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooBernstein.h"

//#include <src/HggOutputReader.h>
#include <src/ZeeOutputReader.h>

#include <iostream>
#include <map>
#include <vector>
using namespace std;
#define nCat 2

#include "selectionMaps.C"


TH2F* getSelectionMap(int map=0,bool isData=true){
  switch(map){
  case 0:
    return getSelectionMap0();
  case 1:
    return getSelectionMap1();
  case 2:
    return getSelectionMap2();
  case 3:
    return getSelectionMap3();
  case 4:
    return getSelectionMap4(isData);
  }
  return 0;
}
int passSelection(TH2F* map,float sigEoE,float etaSC, float pt){
  if(sigEoE < map->GetBinContent(map->FindFixBin(fabs(etaSC),pt))) return 0;
  return 1;
}

bool getBaselineSelection(ZeeOutputReader* h){
  if(h->mass < 80
     || h->mass > 100
     || !h->passloose) return false;
  return true;
}


float rescale(float se, float linLow, float linHigh,float pivot, float offset2){
  if(se < pivot) return linLow*(se-pivot)+offset2;
  else return linHigh*(se-pivot)+offset2;
  return 0;
}

float scaleSEoE(float se,float eta,bool isLead){

  //define the scaling
  // derived Oct 2012
  // 52X MC
  if( fabs(eta) < 1.48 ){
    if(isLead) return rescale(se, 1.047, 1.621, 0.007576, 0.007322);
    else       return rescale(se, 1.101, 1.480, 0.007806, 0.007812);
  }else{
    if(isLead) return rescale(se, 1.136, 1.342, 0.011650, 0.011880);
    else       return rescale(se, 1.590, 1.201, 0.012172, 0.013350);
  }
  return 0;
}

void AddToWorkspace(TString inputFile,RooWorkspace* ws,TString tag, bool isData,int selectionMap=0,bool requireCiC=false,int runMin=0, int runMax=999999){
  //gROOT->ProcessLine(".L /home/amott/HggApp/spin/LeptonSPlots/src/ZeeOutputReader.C+");
  cout << "saveDataSet" <<endl;
  TFile *f = TFile::Open(inputFile);
  TTree *tree = (TTree*)f->Get("HggOutput");
  
  ZeeOutputReader h(tree);
  
  std::cout << "Get Selection Map" << std::endl;
  TH2F* map = getSelectionMap(selectionMap,isData);
  std::cout << map->GetBinContent(1,1) << std::endl;

  // fit variables
  RooRealVar* mass   = new RooRealVar("mass",  "Mass [GeV]", 80., 100.);
  RooRealVar* cosT   = new RooRealVar("cosT",  "cos(theta)", -1, 1);
  RooRealVar* sige1  = new RooRealVar("sigEoE1","#sigma_{E}/E Lead Photon",0,0.1);
  RooRealVar* sige2  = new RooRealVar("sigEoE2","#sigma_{E}/E Lead Photon",0,0.1);
  RooRealVar* evtW   = new RooRealVar("evtWeight","Event Weight",1,0,1e6);

  RooRealVar *totEB = new RooRealVar(Form("%s_EB_totalEvents",tag.Data()),"",0,0,1e9);
  RooRealVar *totEE = new RooRealVar(Form("%s_EE_totalEvents",tag.Data()),"",0,0,1e9);

  std::map<std::pair<int,int>, RooDataSet*> dataMapEB, dataMapEE;
  std::map<std::pair<int,int>, RooDataSet*> *datamap;

  for(int i=0;i<nCat;i++){
    for(int j=0;j<nCat;j++){
      dataMapEB[std::pair<int,int>(i,j)] = new RooDataSet(Form("%s_EB_%d_%d",tag.Data(),i,j),"",RooArgList(*mass,*cosT,*sige1,*sige2,*evtW),"evtWeight");
      dataMapEE[std::pair<int,int>(i,j)] = new RooDataSet(Form("%s_EE_%d_%d",tag.Data(),i,j),"",RooArgList(*mass,*cosT,*sige1,*sige2,*evtW),"evtWeight");
    }
  }

  Long64_t iEntry=-1;
  cout << "Making DataSet" << endl;
  Long64_t nEB=0,nEE=0;
  while(h.GetEntry(++iEntry)){
    if( !(iEntry%10000) ) cout << "Processing " << iEntry <<endl;
    if(!getBaselineSelection(&h)) continue;
    //if(isData && (h.runNumber < runMin || h.runNumber > runMax) ) continue;

    if(fabs(h.Ele1etaSC) < 1.48 && fabs(h.Ele2etaSC) < 1.48) nEB++;
    else nEE++;

    if(requireCiC){
      if(!h.passtight) continue;

    }
    float se1 = h.Ele1sigEscaleoEpho;
    float se2 = h.Ele2sigEscaleoEpho;
    if(!isData){
      //se1 = scaleSEoE(se1,h.Photon_etaSC[1],true);
      //se2 = scaleSEoE(se2,h.Photon_etaSC[0],false);
    }
    int p1 = passSelection(map,se1,h.Ele1etaSC,h.Ele1pt);
    int p2 = passSelection(map,se2,h.Ele2etaSC,h.Ele2pt);
    if(p1 >= nCat || p2 >= nCat) continue;
    datamap = ((fabs(h.Ele1etaSC) < 1.48 && fabs(h.Ele2etaSC) < 1.48) ?
	       &dataMapEB : & dataMapEE);
      
    mass->setVal(h.mass);
    cosT->setVal(0);
    sige1->setVal(se1);
    sige2->setVal(se2);
    evtW->setVal(1);
    (*datamap)[std::pair<int,int>(p1,p2)]->add(RooArgList(*mass,*cosT,*sige1,*sige2));
  }
  cout << "Processed " << iEntry << " Entries" <<endl;

  totEB->setVal(nEB);
  totEE->setVal(nEE);

  std::map<std::pair<int,int>, RooDataSet*>::iterator dIt;
  for(dIt = dataMapEB.begin();dIt!=dataMapEB.end();dIt++) ws->import(*(dIt->second));
  for(dIt = dataMapEE.begin();dIt!=dataMapEE.end();dIt++) ws->import(*(dIt->second));
  ws->import(*totEB);
  ws->import(*totEE);
  cout << "Done" <<endl;
}

void MakeSignalFit(RooWorkspace *ws, TString inputTag, TString outputTag){
  //RooRealVar mass("mass","m_{#gamma#gamma}",90,190);
  RooRealVar mass = *(ws->var("mass"));
  RooRealVar cosT = *(ws->var("cosT"));
  cosT.setBins(5);

  RooRealVar mean(Form("%s_mean",outputTag.Data()),Form("%s_mean",outputTag.Data()),130,100,180);
  RooRealVar sig1(Form("%s_sigma1",outputTag.Data()),Form("%s_sigma1",outputTag.Data()),1,0.2,10);
  RooRealVar sig2(Form("%s_sigma2",outputTag.Data()),Form("%s_sigma2",outputTag.Data()),1,0.2,10);
  RooRealVar f(Form("%s_f",outputTag.Data()),Form("%s_f",outputTag.Data()),0.1,0,1);
  RooGaussian g1(Form("%s_g1",outputTag.Data()),Form("%s_g1",outputTag.Data()),mass,mean,sig1);
  RooGaussian g2(Form("%s_g2",outputTag.Data()),Form("%s_g2",outputTag.Data()),mass,mean,sig2);
  RooAddPdf SignalModel(outputTag.Data(),"Signal Model",RooArgList(g1,g2),f);
  
  RooDataSet *ds = (RooDataSet*)ws->data(inputTag);

  //RooKeysPdf cosTkde(Form("%s_cosTpdf",outputTag.Data()),"KDE for cos(theta) dist",cosT,*ds);
  RooDataHist hist(Form("%s_cosThist",outputTag.Data()),"Data Hist for cos(theta)",RooArgSet(cosT),*ds);
  RooHistPdf cosTkde(Form("%s_cosTpdf",outputTag.Data()),"Hist PDF for cos(theta)",RooArgSet(cosT),hist);

  RooFitResult *res = SignalModel.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2));
  std::cout << res <<std::endl;
  res->SetName(Form("%s_FitResult",outputTag.Data()));

  cosTkde.fitTo(*ds);
  ws->import(SignalModel);
  ws->import(cosTkde);
  ws->import(*res);
}

void MakeSignalFits(RooWorkspace *ws,TString tag){

  for(int iCat1=0;iCat1<nCat;iCat1++){
    for(int iCat2=0;iCat2<nCat;iCat2++){

      MakeSignalFit(ws,Form("%s_EB_%d_%d",tag.Data(),iCat1,iCat2),
		    Form("%s_FIT_EB_%d_%d",tag.Data(),iCat1,iCat2));
      MakeSignalFit(ws,Form("%s_EE_%d_%d",tag.Data(),iCat1,iCat2),
		    Form("%s_FIT_EE_%d_%d",tag.Data(),iCat1,iCat2));

      
    }//iCat2
  }//iCat1

}
void AddSWeight(RooWorkspace *ws, TString tag,TString mcName){

  RooRealVar mean = *(ws->var(Form("%s_FIT_%s_mean",mcName.Data(),tag.Data())));
  mean.setConstant(kTRUE);

  RooDataSet *ds  = (RooDataSet*)ws->data(Form("Data_%s",tag.Data()));
  RooAbsPdf  *pdf = ws->pdf(Form("Data_%s_FIT_%s_fitModel",mcName.Data(),tag.Data()));
  RooRealVar *Ns = ws->var(Form("Data_%s_FIT_%s_Nsig",mcName.Data(),tag.Data()));
  RooRealVar *Nb = ws->var(Form("Data_%s_FIT_%s_Nbkg",mcName.Data(),tag.Data()));

  RooStats::SPlot* sData = new RooStats::SPlot(Form("%s_sData_%s",mcName.Data(),tag.Data()),"",
					       *ds,pdf,RooArgList(*Ns,*Nb));


  RooArgList sweights = sData->GetSWeightVars();
  RooRealVar *Ns_sw = (RooRealVar*)sweights.at(0); //ws->var(Form("Data_%s_FIT_%s_Nsig_sw",mcName.Data(),tag.Data()));
  RooRealVar *Nb_sw = (RooRealVar*)sweights.at(1); //ws->var(Form("Data_%s_FIT_%s_Nbkg_sw",mcName.Data(),tag.Data()));

  //make weighted datasets
  RooRealVar *cosT  = ws->var("cosT");
  RooRealVar *sige1 = ws->var("sigEoE1");
  RooRealVar *sige2 = ws->var("sigEoE2");
  RooDataSet sig(Form("%s_sigWeight_%s",mcName.Data(),tag.Data()),"",ds,RooArgSet(*cosT,*sige1,*sige2,*Ns_sw,*Nb_sw),0,Ns_sw->GetName());
  RooDataSet bkg(Form("%s_bkgWeight_%s",mcName.Data(),tag.Data()),"",ds,RooArgSet(*cosT,*sige1,*sige2,*Ns_sw,*Nb_sw),0,Nb_sw->GetName());
  

  ws->import(*ds,RooFit::Rename(Form("%s_sData_%s",mcName.Data(),tag.Data())));
  ws->import(sig);
  ws->import(bkg);

}



void MakeBackgroundFit(RooWorkspace *ws,TString tag, TString mcName,float initMass,float range,bool gausPen){
  RooRealVar mass = *(ws->var("mass"));
  cout << "1" <<endl;
  //build the signal model
  float min = initMass - (gausPen ? 3 : 1)*range; //allow 3 sigma up down
  float max = initMass + (gausPen ? 3 : 1)*range; //allow 3 sigma up down
  RooRealVar mean(Form("Data_%s_FIT_%s_mean",mcName.Data(), tag.Data()),"mean",initMass,min,max);
  cout << "5" <<endl;
  //won't use this unless gaussian penalty is on
  RooRealVar initialMass(Form("Data_%s_FIT_%s_initialMass",mcName.Data(), tag.Data()),"Initial mass val",
			 initMass);
  cout << "3" <<endl;
  RooRealVar penaltyMass(Form("Data_%s_FIT_%s_penaltyMass",mcName.Data(), tag.Data()),"Mass Penalty",
			 range);
  cout << "4" <<endl;
  RooGaussian meanConst("meanConstraint","",mean,initialMass,penaltyMass);
  cout << "5   " << Form("%s_FIT_%s_sigma1",mcName.Data(),tag.Data()) << endl;
  RooRealVar sig1Base = *(ws->var(Form("%s_FIT_%s_sigma1",mcName.Data(),tag.Data())));
  cout << "6" <<endl;
  RooRealVar sig2Base = *(ws->var(Form("%s_FIT_%s_sigma2",mcName.Data(),tag.Data())));
  cout << "7" <<endl;
  RooRealVar fBase    = *(ws->var(Form("%s_FIT_%s_f",mcName.Data(),tag.Data())));


  RooRealVar sig1Width("sig1Wdith","",sig1Base.getError()/8);
  RooRealVar sig2Width("sig2Wdith","",sig2Base.getError()/8);
  RooRealVar fWidth("fWdith","",fBase.getError()/8);

  RooRealVar sig1(Form("Data_%s_FIT_%s_sigma1",mcName.Data(), tag.Data()),"sig1",sig1Base.getVal(),sig1Base.getVal()-4*sig1Base.getError(),sig1Base.getVal()+4*sig1Base.getError());
  RooRealVar sig2(Form("Data_%s_FIT_%s_sigma2",mcName.Data(), tag.Data()),"sig2",sig2Base.getVal(),sig2Base.getVal()-4*sig2Base.getError(),sig2Base.getVal()+4*sig2Base.getError());
  RooRealVar f(Form("Data_%s_FIT_%s_f",mcName.Data(), tag.Data()),"f",fBase.getVal(),fBase.getVal()-4*fBase.getError(),fBase.getVal()+4*fBase.getError());
  
  RooGaussian sig1Const("sig1Constraint","",sig1,sig1Base,sig1Width);
  RooGaussian sig2Const("sig2Constraint","",sig2,sig2Base,sig2Width);
  RooGaussian fConst("fConstraint","",f,fBase,fWidth);
  //sig1.setConstant(kTRUE);
  //sig2.setConstant(kTRUE);
  //f.setConstant(kTRUE);
  RooGaussian g1(Form("Data_%s_FIT_%s_g1",mcName.Data(),tag.Data()),"g1",mass,mean,sig1);
  RooGaussian g2(Form("Data_%s_FIT_%s_g2",mcName.Data(),tag.Data()),"g2",mass,mean,sig2);
  RooAddPdf SignalModel(Form("Data_%s_FIT_%s_signalModel",mcName.Data(),tag.Data()),
			"Signal Model",RooArgList(g1,g2),f);

RooAbsData *ds = ws->data(Form("Data_%s",tag.Data()));

  //background model
  cout << "8" <<endl;

  
  RooRealVar alpha1(Form("Data_%s_FIT_%s_alpha1",mcName.Data(),tag.Data()),"alpha1",-0.1,-1.,0.);
  RooRealVar alpha2(Form("Data_%s_FIT_%s_alpha2",mcName.Data(),tag.Data()),"alpha2",-0.1,-1.,0.);
  RooRealVar f_bkg( Form("Data_%s_FIT_%s_f",mcName.Data(),tag.Data()),"f_bkg",0.1,0,1);
  RooExponential exp1(Form("Data_%s_FIT_%s_exp1",mcName.Data(),tag.Data()),"exp1",mass,alpha1);
  RooExponential exp2(Form("Data_%s_FIT_%s_exp2",mcName.Data(),tag.Data()),"exp2",mass,alpha2);
    
  RooAddPdf BkgModel(Form("Data_%s_FIT_%s_bkgModel",mcName.Data(),tag.Data()),"Background Model",
		     RooArgList(exp1,exp2),f_bkg);
  
  /*
  RooRealVar pC(Form("Data_%s_FIT_%s_pC",mcName.Data(),tag.Data()),"pC",1,1,1);
  RooRealVar p0(Form("Data_%s_FIT_%s_p0",mcName.Data(),tag.Data()),"p0",0,-1e6,1e6);
  RooRealVar p1(Form("Data_%s_FIT_%s_p1",mcName.Data(),tag.Data()),"p1",0,-1e6,1e6);
  RooRealVar p2(Form("Data_%s_FIT_%s_p2",mcName.Data(),tag.Data()),"p2",0,-1e6,1e6);
  RooRealVar p3(Form("Data_%s_FIT_%s_p3",mcName.Data(),tag.Data()),"p3",0,-1e6,1e6);

  RooFormulaVar pCmod(Form("Data_%s_FIT_%s_pCmod",mcName.Data(),tag.Data()),"","@0*@0",pC);
  RooFormulaVar p0mod(Form("Data_%s_FIT_%s_p0mod",mcName.Data(),tag.Data()),"","@0*@0",p0);
  RooFormulaVar p1mod(Form("Data_%s_FIT_%s_p1mod",mcName.Data(),tag.Data()),"","@0*@0",p1);
  RooFormulaVar p2mod(Form("Data_%s_FIT_%s_p2mod",mcName.Data(),tag.Data()),"","@0*@0",p2);
  RooFormulaVar p3mod(Form("Data_%s_FIT_%s_p3mod",mcName.Data(),tag.Data()),"","@0*@0",p3);

  

  RooArgList *args;
  if(ds->sumEntries() < 300) args = new RooArgList(pCmod,p0mod,p1mod,p2mod);
  else args = new RooArgList(pCmod,p0mod,p1mod,p2mod,p3mod);

  RooBernstein BkgModel(Form("Data_%s_FIT_%s_bkgModel",mcName.Data(),tag.Data()),"Background Model",mass,*args);
  */
  RooRealVar Nsig(Form("Data_%s_FIT_%s_Nsig",mcName.Data(),tag.Data()),"N signal Events",100,0,100000);
  RooRealVar Nbkg(Form("Data_%s_FIT_%s_Nbkg",mcName.Data(),tag.Data()),"N background Events",100,0,1e+09);

  //fit model
  RooAddPdf FitModel(Form("Data_%s_FIT_%s_fitModel",mcName.Data(),tag.Data()),"Fit Model",
		     RooArgList(SignalModel,BkgModel),RooArgList(Nsig,Nbkg));

  //do the fit

  RooArgSet allConstraints(sig1Const,sig2Const,fConst);
  if(gausPen) allConstraints.add(meanConst);
  


  FitModel.fitTo(*ds,RooFit::ExternalConstraints(allConstraints),RooFit::Strategy(0),RooFit::Minos(kFALSE));
  RooFitResult *res=FitModel.fitTo(*ds,RooFit::ExternalConstraints(allConstraints),RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::Minos(kFALSE));
  res->SetName(Form("Data_%s_FIT_%s_fitResult",mcName.Data(),tag.Data()) );

  ws->import(FitModel);
  ws->import(*res);

  AddSWeight(ws,tag,mcName);
}

void makeBackgroundFits(RooWorkspace *ws){
  MakeBackgroundFit(ws, "EB_0_0","DY",124.5,2,false); // do the initial fit
  float massPeak  = ws->var("Data_DY_FIT_EB_0_0_mean")->getVal();
  cout << massPeak << endl;
  float massRange = ws->var("Data_DY_FIT_EB_0_0_mean")->getError();
  cout << massRange <<endl;
  for(int iCat1=0;iCat1<nCat;iCat1++){
    for(int iCat2=0;iCat2<nCat;iCat2++){
      MakeBackgroundFit(ws, Form("EE_%d_%d",iCat1,iCat2),"DY",massPeak,massRange,true); // do the initial fit      
      if(iCat1==0 && iCat2==0) continue;
      MakeBackgroundFit(ws, Form("EB_%d_%d",iCat1,iCat2),"DY",massPeak,massRange,true); // do the initial fit      
      
    }//iCat2
  }//iCat1

}

void MakeWorkspaces(TString inputDataFile,TString inputHggFile,
		    TString outputFile,int selectionMap=0,bool requireCiC=false,int runMin=0,int runMax=999999){
  
  TFile *out = TFile::Open(outputFile,"RECREATE");
  RooWorkspace * ws = new RooWorkspace();
  ws->SetName("cms_zee_spin_workspace");  

  AddToWorkspace(inputDataFile,ws,"Data",true,selectionMap,requireCiC,runMin,runMax);
  AddToWorkspace(inputHggFile, ws,"DY",false,selectionMap,requireCiC);
  
  out->cd();
  ws->Write();
  out->Close();
}

void RunAllFits(TString dataWorkspaceFile,TString fitWorkspaceFile){
  TFile *f = TFile::Open(dataWorkspaceFile);
  RooWorkspace *ws = (RooWorkspace*)f->Get("cms_zee_spin_workspace");
  MakeSignalFits(ws,"DY");
  makeBackgroundFits(ws);
  
  TFile *out = TFile::Open(fitWorkspaceFile,"RECREATE");
  out->cd();
  ws->Write();
  out->Close();
  f->Close();
}

TCanvas* DrawFit(TString wsFile, TString tag, TString mcType){
  TCanvas *cv = new TCanvas(Form("%s_%s",mcType.Data(),tag.Data()));

  TFile *f = TFile::Open(wsFile);
  RooWorkspace *ws = (RooWorkspace*)f->Get("cms_zee_spin_workspace");
  
  RooRealVar* mass = ws->var("mass");
  RooPlot* frame  = mass->frame(100,180,80);

  ws->data(Form("Data_%s",tag.Data()))->plotOn(frame);

  double Ns = ws->var(Form("Data_%s_FIT_%s_Nsig",mcType.Data(),tag.Data()))->getVal();
  double Nb = ws->var(Form("Data_%s_FIT_%s_Nbkg",mcType.Data(),tag.Data()))->getVal();
  
  ws->pdf(Form("Data_%s_FIT_%s_fitModel",mcType.Data(),tag.Data()))->plotOn(frame, RooFit::LineColor(kBlack));
  ws->pdf(Form("Data_%s_FIT_%s_bkgModel",mcType.Data(),tag.Data()))->plotOn(frame, RooFit::Normalization(Nb/(Ns+Nb)),RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));
  
  frame->Draw();
  return cv;
}



void ListCats(TString wsFile,TString mcType){
  TFile *f = TFile::Open(wsFile);
  RooWorkspace *ws = (RooWorkspace*)f->Get("cms_zee_spin_workspace");
 
  cout << "\n**\t\tEB\t\t**\n"<<endl;
  for(int i=0;i<nCat;i++){
    for(int j=0;j<nCat;j++){
      double Ns  = ws->var(Form("Data_%s_FIT_EB_%d_%d_Nsig",mcType.Data(),i,j))->getVal();
      double Nse = ws->var(Form("Data_%s_FIT_EB_%d_%d_Nsig",mcType.Data(),i,j))->getError();
      if(Ns > Nse/3) cout << i << "  " << j << "\t" << Ns << " +/- " << Nse <<endl; 
    }
  }

  cout << "\n**\t\tEE\t\t**\n"<<endl;
  for(int i=0;i<nCat;i++){
    for(int j=0;j<nCat;j++){
      double Ns  = ws->var(Form("Data_%s_FIT_EE_%d_%d_Nsig",mcType.Data(),i,j))->getVal();
      double Nse = ws->var(Form("Data_%s_FIT_EE_%d_%d_Nsig",mcType.Data(),i,j))->getError();
      if(Ns > Nse/3) cout << i << "  " << j << "\t" << Ns << " +/- " << Nse <<endl; 
    }
  }
}

