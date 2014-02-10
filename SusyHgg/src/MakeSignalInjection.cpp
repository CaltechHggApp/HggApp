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

#include "TFile.h"

#include "../include/ArgParser.hh"

#define NUM_CPU 1

RooAbsPdf* getBackgroundFit(RooAbsData& data,RooRealVar& mgg);
RooDataSet* getToyDataSet(RooAbsData* data, RooRealVar& mgg,bool isData,float scale_factor=0);

int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("InputDataFile","Input file containing the data");
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
  TFile inputSMSFile(a.getArgument("InputSMSFile").c_str());

  RooWorkspace* inputDataWs = (RooWorkspace*)inputDataFile.Get("susy_hgg_workspace");
  RooWorkspace* inputSMSWs  = (RooWorkspace*)inputSMSFile.Get("susy_hgg_workspace");


  RooWorkspace* outputWs = new RooWorkspace("susy_hgg_workspace","susy_hgg_workspace");
  
  //import SM Higgs datasets
  outputWs->import( *inputDataWs->data("gg_H_125") );
  outputWs->import( *inputDataWs->data("vbf_H_125") );
  outputWs->import( *inputDataWs->data("wz_H_125") );
  
  RooDataSet * inputData = (RooDataSet*)inputDataWs->data("data");
  RooRealVar * mgg       = inputDataWs->var("mgg");

  RooDataSet* toyData = getToyDataSet(inputData,*mgg,true);
  toyData->SetName("data");

  float lumi = 19780;
  if(a.longFlagPres("lumi")) lumi = atof(a.getLongFlag("lumi").c_str());
  float sms_xsec = 1;
  if(a.longFlagPres("crossSection")) sms_xsec = atof(a.getLongFlag("crossSection").c_str());
  


  toyData->append( *getToyDataSet(inputDataWs->data("gg_H_125"),*mgg,false, 19.27*2.28E-3*lumi/inputDataWs->var("N_gg_H_125")->getVal()) );
  toyData->append( *getToyDataSet(inputDataWs->data("vbf_H_125"),*mgg,false, 1.578*2.28E-3*lumi/inputDataWs->var("N_vbf_H_125")->getVal()) );
  toyData->append( *getToyDataSet(inputDataWs->data("wz_H_125"),*mgg,false, 0.7046*2.28E-3*lumi/inputDataWs->var("N_wz_H_125")->getVal()) );

  if(sms_xsec*2.28E-3*lumi/inputSMSWs->var(Form("N_%s",a.getArgument("SMSName").c_str()))->getVal() >=0.6)
    toyData->append( *getToyDataSet(inputSMSWs->data(a.getArgument("SMSName").c_str()),*mgg,false,sms_xsec*2.28E-3*lumi/inputSMSWs->var(Form("N_%s",a.getArgument("SMSName").c_str()))->getVal()) );

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


RooAbsPdf* getBackgroundFit(RooAbsData& data,RooRealVar& mgg) {
  RooRealVar* a1 = new RooRealVar("alpha1","",-0.1,-1.,0.);
  RooRealVar* a2 = new RooRealVar("alpha2","",-0.1,-1.,0.);
  RooRealVar* f  = new RooRealVar("f","",0.1,0.,1.);

  RooExponential* exp1 = new RooExponential("exp1","",mgg,*a1);
  RooExponential* exp2 = new RooExponential("exp2","",mgg,*a2);
  RooAddPdf     * dexp = new RooAddPdf("dexp","",RooArgList(*exp1,*exp2),*f);

  RooRealVar *Nbkg = new RooRealVar("Nbkg","",1000,0,1e9);

  RooExtendPdf *bkgModel = new RooExtendPdf("bkgModel","",*dexp,*Nbkg);

  
  bkgModel->fitTo(data,RooFit::Save(kFALSE),RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU));
  bkgModel->fitTo(data,RooFit::Save(kFALSE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU));

  return bkgModel;
}

RooDataSet* getToyDataSet(RooAbsData* data, RooRealVar& mgg,bool isData,float scale_factor) {
  RooAbsPdf* mggPdf=0;
  if(isData) {
    mggPdf = getBackgroundFit(*data,mgg);
  }else{ //its an SMS, so use a finally binned histogram
    mgg.setBins(80);
    RooDataHist *hist = new RooDataHist("hist","",mgg,*data);
    mggPdf = new RooHistPdf("mggPdf","",mgg,*hist);
  }

  int N = (isData ? data->sumEntries() : data->sumEntries()*scale_factor);

  std::cout << ">>>>>>>>>>>> N=" << N <<  std::endl;

  RooDataSet* toyData =  mggPdf->generate(mgg, N, RooFit::Extended(), RooFit::ProtoData(*(RooDataSet*)data,true));  
  
  delete mggPdf;
  return toyData;
}
