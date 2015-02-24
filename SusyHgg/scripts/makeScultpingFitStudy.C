#include "makeMggFit.C"


#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

#include <cfloat>
#include <iostream>
#include "assert.h"

RooWorkspace* makeFit(TTree* tree,TString wsName,float mggMax=160);

void computeSF(RooWorkspace* w, TString pdfName, TString sfName, float sigEff);

void makeSculptingFitStudy(TString outputFileName) {

  const int nBins=3;
  TString labels[nBins] = {"HighPt","HighRes","LowRes"};
  TString dataFiles[nBins] = {
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v8__HighPt.root",
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v8__HighRes.root",
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v8__LowRes.root"
  };

  float sigEffs[nBins] = {1.52,1.48,2.5};

  TFile outputFile(outputFileName,"RECREATE");

  

  for(int iBox=0;iBox<nBins;iBox++) {
    TFile inputFile(dataFiles[iBox]);
    TTree *t = (TTree*)inputFile.Get("SusyHggTree");
    for(int end=150;end<=250;end+=10) {
      RooWorkspace *wa = makeFit(t,Form("%s_ws_allMR_%d",labels[iBox].Data(),end),end);
      RooWorkspace *wl = makeFit(t->CopyTree("MR>=250 && MR<400"),Form("%s_ws_lowMR_%d",labels[iBox].Data(),end),end);
      RooWorkspace *wh = makeFit(t->CopyTree("MR>=800"),Form("%s_ws_highMR_%d",labels[iBox].Data(),end),end);

      outputFile.cd();

      computeSF(wa,"dexp","sf_dexp",sigEffs[iBox]);
      computeSF(wl,"dexp","sf_dexp",sigEffs[iBox]);
      computeSF(wh,"dexp","sf_dexp",sigEffs[iBox]);

      wa->Write();
      wl->Write();
      wh->Write();
  }
 
  }

}



RooWorkspace* makeFit(TTree* tree,TString wsName, float mggMax) {
  RooWorkspace *ws = new RooWorkspace(wsName,"");

  std::vector< TString (*)(TString, RooRealVar&, RooWorkspace&) > bkgPdfList;
  bkgPdfList.push_back(makeDoubleExp);
 //  bkgPdfList.push_back(makeSingleExp);
//   bkgPdfList.push_back(makeModExp);
//   bkgPdfList.push_back(makeSinglePow);
//   bkgPdfList.push_back(makeDoublePow);



  RooRealVar mgg("mgg","m_{#gamma#gamma}",103,mggMax,"GeV");
  mgg.setBins(mggMax-103);

  mgg.setRange("sideband_low", 103,120);
  mgg.setRange("sideband_high",131,160);

  RooRealVar MR("MR","",0,3000,"GeV");
  MR.setBins(60);
  
  RooRealVar Rsq("Rsq","",0,1,"GeV");
  Rsq.setBins(20);

  RooRealVar hem1_M("hem1_M","",-1,2000,"GeV");
  hem1_M.setBins(40);

  RooRealVar hem2_M("hem2_M","",-1,2000,"GeV");
  hem2_M.setBins(40);

  RooRealVar ptgg("ptgg","p_{T}^{#gamma#gamma}",0,500,"GeV");
  ptgg.setBins(50);

  RooDataSet data("data","",tree,RooArgSet(mgg,MR,Rsq,hem1_M,hem2_M,ptgg));
  for(auto func = bkgPdfList.begin(); func != bkgPdfList.end(); func++) {
    TString tag = (*func)("bkgPdf",mgg,*ws);
    ws->pdf("bkgPdf_"+tag+"_ext")->fitTo(data,RooFit::Strategy(0),RooFit::Extended(kTRUE));
    RooFitResult* bres = ws->pdf("bkgPdf_"+tag+"_ext")->fitTo(data,RooFit::Strategy(2),RooFit::Save(kTRUE),RooFit::Extended(kTRUE));
    bres->Print("V");
    assert(bres->covQual()>=2);

    
    bres->SetName("fitres_bkgPdf_"+tag);
    ws->import(*bres);
  }

  ws->import(data);

  return ws;
}

void computeSF(RooWorkspace* w, TString pdfName, TString sfName, float sigEff) {
  RooAbsPdf * pdf = w->pdf("bkgPdf_"+pdfName+"_ext");
  RooRealVar* mgg = w->var("mgg");
  RooFitResult* res = (RooFitResult*)w->obj("fitres_bkgPdf_"+pdfName);

  mgg->setRange("sig",125-2*sigEff, 126+2*sigEff);

  RooRealVar* Nbkg = w->var("bkgPdf_dexp_Nbkg");

  float nLow  = w->data("data")->sumEntries("mgg>=103 && mgg<120");
  float nHigh = w->data("data")->sumEntries("mgg>=131 && mgg<160");

  RooAbsReal* integral = pdf->createIntegral(*mgg,RooFit::NormSet(*mgg),RooFit::Range("sig"));
  std::cout << integral->getVal() << "  " << integral->getPropagatedError(*res) << std::endl;
  float sf = integral->getVal()*Nbkg->getVal()/(nLow+nHigh);

  float sfe = sf*sqrt( 1/(nLow+nHigh)+pow( integral->getPropagatedError(*res)/integral->getVal(),2 ) + 1/Nbkg->getVal() );

  RooRealVar sfRV(sfName,"",sf);
  sfRV.setError(sfe);

  w->import(sfRV);
}
