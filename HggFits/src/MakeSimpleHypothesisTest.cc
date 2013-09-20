#include "MakeSimpleHypothesisTest.h"

#include "assert.h"
#include <iostream>

MakeSimpleHypothesisTest::MakeSimpleHypothesisTest(const TString& inputFileName, const TString& outputFileName) {
  if(inputFileName != ""){
    //opens the input file and gets the input workspace
      inputFile = new TFile(inputFileName);
      ws = ((RooWorkspace*)inputFile->Get("cms_hgg_spin_workspace"));
      //ws->Print();
      //extract MC labels from the input workspace
      //MakeSpinFits::getLabels("labels",&mcLabel,ws);
      //extract category labels from the input workspace
      MakeSpinFits::getLabels("evtcat",&catLabels,ws);
  }
  if(outputFileName != ""){
      //opens the output file
      outputFile = new TFile(outputFileName,"RECREATE");
    outputFile->cd();
    ws->Write(ws->GetName(),TObject::kWriteDelete);
  }
}

void MakeSimpleHypothesisTest::MakeFakeSignalModels(const TString& baseMCName,const TString& tag,float scaleFactor,float  meanVal) {
    assert(ws!=0);
    std::cout << Form("%s_FIT_%s_sigma1",baseMCName.Data(),tag.Data()) << std::endl;
    RooRealVar* mass = ws->var("mass");

    RooRealVar mean( Form("Data_Hgg%.0f_FIT_%s_mean",meanVal,tag.Data()),"", meanVal);
    RooRealVar sig1( Form("Data_Hgg%.0f_FIT_%s_sigma1",meanVal,tag.Data()),"",
                     ws->var( Form("%s_FIT_%s_sigma1",baseMCName.Data(),tag.Data()) )->getVal()*scaleFactor );
    RooRealVar sig2( Form("Data_Hgg%.0f_FIT_%s_sigma2",meanVal,tag.Data()),"",
                     ws->var( Form("%s_FIT_%s_sigma2",baseMCName.Data(),tag.Data()) )->getVal()*scaleFactor );
    RooRealVar sig3( Form("Data_Hgg%.0f_FIT_%s_sigma3",meanVal,tag.Data()),"",
                     ws->var( Form("%s_FIT_%s_sigma3",baseMCName.Data(),tag.Data()) )->getVal()*scaleFactor );


    RooRealVar f1( Form("Data_Hgg%.0f_FIT_%s_f1",meanVal,tag.Data()),"",
                   ws->var( Form("%s_FIT_%s_f1",baseMCName.Data(),tag.Data()) )->getVal() );
    RooRealVar f2( Form("Data_Hgg%.0f_FIT_%s_f2",meanVal,tag.Data()),"",
                   ws->var( Form("%s_FIT_%s_f2",baseMCName.Data(),tag.Data()) )->getVal() );

    mean.setConstant(kTRUE);
    sig1.setConstant(kTRUE);
    sig2.setConstant(kTRUE);
    sig2.setConstant(kTRUE);
    f1.setConstant(kTRUE);
    f2.setConstant(kTRUE);

    RooGaussian g1( Form("Data_Hgg%.0f_FIT_%s_g1",meanVal,tag.Data()), "",*mass,mean,sig1); 
    RooGaussian g2( Form("Data_Hgg%.0f_FIT_%s_g2",meanVal,tag.Data()), "",*mass,mean,sig2); 
    RooGaussian g3( Form("Data_Hgg%.0f_FIT_%s_g3",meanVal,tag.Data()), "",*mass,mean,sig3); 

    RooAddPdf SignalModel (Form("Data_Hgg%.0f_FIT_%s",meanVal,tag.Data()),"Signal Model",RooArgList(g1,g2,g3),RooArgList(f1,f2));

    ws->import(SignalModel);
}


void MakeSimpleHypothesisTest::MakeAllFakeSignalModels() {
    for(int mh=mh_low; mh<mh_high; mh+=mh_step) { //loop over the masses        
        for(auto catIt = catLabels.begin(); catIt !=catLabels.end(); catIt++) { //loop over the categories
            std::cout << *catIt << std::endl;
            MakeFakeSignalModels("jhu0plus125",*catIt,mh/125.,float(mh));
        }        
    }
}

void MakeSimpleHypothesisTest::MakeHypothesisTest() {
    TH1F significance("CombinedSignificance","",401,100,500);
    TH1F NSig("signalYield","",401,100,500);

    MakeSpinFits fitter("","");
    fitter.setWorkspace(ws);
    fitter.setEmulatedMassHack("jhu0plus125");
    
    for(int mh=mh_low; mh<mh_step; mh+=mh_step) { //loop over the masses
        fitter.MakeCombinedSignalTest( Form("Hgg%d",mh),false);
        float NS  = ws->var( Form("Data_Hgg%d_FULLFIT_Nsig",mh) )->getVal();
        float NSE = ws->var( Form("Data_Hgg%d_FULLFIT_Nsig",mh) )->getError();
        NSig.SetBinContent(mh-125,NS);
        significance.SetBinContent(mh-125,NS/NSE);
    }
    outputFile->cd();
    NSig.Write();
    significance.Write();
    outputFile->Close();
    
}

void MakeSimpleHypothesisTest::run() {
    MakeAllFakeSignalModels();
    MakeHypothesisTest();
}
