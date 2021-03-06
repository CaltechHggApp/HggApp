#include <iostream>
#include <string>

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooCategory.h"

#include "TFile.h"

#include "../include/ArgParser.hh"

#define NUM_CPU 1

RooDataSet* getToyDataSet(RooAbsData* data,RooRealVar& mgg,bool isData,RooAbsPdf* genPdf=0,float scale_factor=0);

int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("InputDataFile","Input file containing the data");
  a.addArgument("InputDataFitFile","Input file containing the data fits");
  a.addArgument("InputSMSFile","Input file containing the SMS");
  a.addArgument("SMSName","Name of the SMS to inject");
  a.addArgument("OutputFileName","Name of the Output File");

  a.addLongOption("crossSection",ArgParser::reqArg,"x-section for the SMS (in pb) [def: 1pb]");

  a.addLongOption("lumi",ArgParser::reqArg,"data luminosity (in pb^-1) [def: 19780 pb^-1]");

  std::string ret;
  if(a.process(ret) !=0){
    std::cout << "Invalid Options:  " << ret << std::endl;
    a.printOptions(argv[0]);
    return 0;
  }

  TFile inputDataFile(a.getArgument("InputDataFile").c_str());
  TFile inputDataFitFile(a.getArgument("InputDataFitFile").c_str());
  TFile inputSMSFile(a.getArgument("InputSMSFile").c_str());

  RooWorkspace* inputDataWs = (RooWorkspace*)inputDataFile.Get("susy_hgg_workspace");
  RooWorkspace* inputDataFitWs = (RooWorkspace*)inputDataFitFile.Get("susy_hgg_workspace");
  RooWorkspace* inputSMSWs  = (RooWorkspace*)inputSMSFile.Get("susy_hgg_workspace");

  RooCategory* cats = (RooCategory*)inputDataFitWs->obj("evtcat");

  RooWorkspace* outputWs = new RooWorkspace("susy_hgg_workspace","susy_hgg_workspace");
  
  //import SM Higgs datasets
  outputWs->import( *inputDataWs->data("gg_H_125") );
  outputWs->import( *inputDataWs->data("vbf_H_125") );
  outputWs->import( *inputDataWs->data("wz_H_125") );
  
  RooDataSet * inputData = (RooDataSet*)inputDataFitWs->data("data_Combined");
  RooRealVar * mgg       = inputDataWs->var("mgg");

  float lumi = 19780;
  if(a.longFlagPres("lumi")) lumi = atof(a.getLongFlag("lumi").c_str());
  float sms_xsec = 1;
  if(a.longFlagPres("crossSection")) sms_xsec = atof(a.getLongFlag("crossSection").c_str());
  
  RooDataSet* toyData=0;

  int catIndex=-1;
  while(cats->setIndex(++catIndex)==0) {

    TString catName = cats->getLabel();

    RooAbsPdf* bkgPdf = inputDataFitWs->pdf(catName+"_bkgModel");

    if(toyData==0) {
      toyData = getToyDataSet(inputData->reduce("evtcat==evtcat::"+catName),*mgg,true,bkgPdf);
      toyData->SetName("data");
    }else{
      toyData->append( *getToyDataSet(inputData->reduce("evtcat==evtcat::"+catName),*mgg,true,bkgPdf) );
    }

//     toyData->append( *getToyDataSet(inputDataFitWs->data("gg_H_125_Combined")->reduce("evtcat==evtcat::"+catName),*mgg,false,0, 19.27*2.28E-3*lumi/inputDataWs->var("N_gg_H_125")->getVal()) );
//     toyData->append( *getToyDataSet(inputDataFitWs->data("vbf_H_125_Combined")->reduce("evtcat==evtcat::"+catName),*mgg,false,0, 1.578*2.28E-3*lumi/inputDataWs->var("N_vbf_H_125")->getVal()) );
//     toyData->append( *getToyDataSet(inputDataFitWs->data("wz_H_125_Combined")->reduce("evtcat==evtcat::"+catName),*mgg,false,0, 0.7046*2.28E-3*lumi/inputDataWs->var("N_wz_H_125")->getVal()) );

    toyData->Print("V");
    RooDataSet *toySMHiggs = getToyDataSet(inputDataFitWs->data("SMTot_Combined")->reduce("evtcat==evtcat::"+catName),*mgg,false,0,1);
    toyData->append( *toySMHiggs );
    //if(sms_xsec*2.28E-3*lumi/inputSMSWs->var(Form("N_%s",a.getArgument("SMSName").c_str()))->getVal() >=0.6) {
    if( inputSMSWs->data(Form("sms_%s_Combined",a.getArgument("SMSName").c_str()))->sumEntries()*sms_xsec >=0.6 ){
      std::cout << "INCLUDING SMS" << std::endl;
      toyData->append( *getToyDataSet(inputSMSWs->data(Form("sms_%s_Combined",a.getArgument("SMSName").c_str()))->reduce("evtcat==evtcat::"+catName),*mgg,false,0,sms_xsec*2.28E-3*lumi/inputSMSWs->var(Form("N_sms_%s",a.getArgument("SMSName").c_str()))->getVal()) );
      //toyData->append( *getToyDataSet(inputSMSWs->data(Form("sms_%s_Combined",a.getArgument("SMSName").c_str()))->reduce("evtcat==evtcat::"+catName),*mgg,false,0,sms_xsec) );
    }else{
      std::cout << "EXCLUDING SMS (x-sec too small)!" << std::endl;
    }
  }
  outputWs->import(*toyData);

  outputWs->import(*inputDataWs->var("N_data"));
  outputWs->import(*inputDataWs->var("N_gg_H_125"));
  outputWs->import(*inputDataWs->var("N_vbf_H_125"));
  outputWs->import(*inputDataWs->var("N_wz_H_125"));
  
  TFile outputFile(a.getArgument("OutputFileName").c_str(),"RECREATE");
  outputWs->Write();

  inputDataFile.Get("mc_name_array")->Write("mc_name_array",TObject::kSingleKey);
  outputFile.Close();

  return 0;
}

RooDataSet* getToyDataSet(RooAbsData* data, RooRealVar& mgg,bool isData,RooAbsPdf* genPdf,float scale_factor) {
  std::cout << "getToyDataSet" << std::endl;
  if(genPdf==0) { //if we don't supply a generator pdf, use a binned one
    mgg.setBins(80);
    RooDataHist *hist = new RooDataHist("hist","",mgg,*data);
    genPdf = new RooHistPdf("mggPdf","",mgg,*hist);
  }

  int N = (isData ? data->sumEntries() : data->sumEntries()*scale_factor);

  std::cout << ">>>>>>>>>>>> N=" << N <<  std::endl;

  //RooDataSet* toyData =  mggPdf->generate(mgg, N, RooFit::Extended(), RooFit::ProtoData(*(RooDataSet*)data,true));  
  RooDataSet* toyData =  genPdf->generate(mgg, RooFit::NumEvents(N+1), RooFit::ProtoData(*(RooDataSet*)data,true,true));  
  
  assert(toyData!=0);

  toyData->Print("V");

  return toyData;
}
