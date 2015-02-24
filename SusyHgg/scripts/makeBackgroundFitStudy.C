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

TGraphAsymmErrors* fitChoiceSystematic(RooWorkspace* w, const int nPdfs, const TString * pdfNames, const int *nPars, TString nominal="dexp");

void computeSigInt(RooWorkspace* w, TString pdfName, TString sfName, float sigEff);

void makeBackgroundFitStudy(TString outputFileName) {

  const int nBins=5;
  TString labels[nBins] = {"HighPt","Hbb","HighRes","LowRes","Zbb"};
  TString dataFiles[nBins] = {
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v7__HighPt.root",
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v7__Hbb.root",
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v7__HighRes.root",
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v7__LowRes.root",
    "DATA/DoublePhoton_22Jan2013_Run2012ABCD_v7__Zbb.root"
  };

  float sigEffs[nBins] = {1.52,2,2,1.48,2.5};

  const int nPdfs = 5;
  TString pdfNames[nPdfs] = {"dexp","sexp","mexp","spow","dpow"};
  int nPar[nPdfs] = {3,1,2,1,3};
  TFile outputFile(outputFileName,"RECREATE");

  

  for(int iBin=0;iBin<nBins;iBin++) {
    TFile inputFile(dataFiles[iBin]);
    TTree *t = (TTree*)inputFile.Get("SusyHggTree");
    RooWorkspace *w160 = makeFit(t,labels[iBin]+"_ws160");
    RooWorkspace *w180 = makeFit(t,labels[iBin]+"_ws180",180);

    outputFile.cd();

    for(int iPdf=0;iPdf<nPdfs;iPdf++) {
      computeSigInt(w160,pdfNames[iPdf],"sigInt_160_"+pdfNames[iPdf],sigEffs[iBin]);
      computeSigInt(w180,pdfNames[iPdf],"sigInt_180_"+pdfNames[iPdf],sigEffs[iBin]);
    }

    w160->Write();
    w180->Write();

    TGraphAsymmErrors* e160 = fitChoiceSystematic(w160,nPdfs,pdfNames,nPar);
    TGraphAsymmErrors* e180 = fitChoiceSystematic(w180,nPdfs,pdfNames,nPar);

    e160->SetName("FitChoice_syst_160_"+labels[iBin]);
    e180->SetName("FitChoice_syst_180_"+labels[iBin]);

    e160->Write();
    e180->Write();
  }
 
};



RooWorkspace* makeFit(TTree* tree,TString wsName, float mggMax) {
  RooWorkspace *ws = new RooWorkspace(wsName,"");

  std::vector< TString (*)(TString, RooRealVar&, RooWorkspace&) > bkgPdfList;
  bkgPdfList.push_back(makeDoubleExp);
  bkgPdfList.push_back(makeSingleExp);
  bkgPdfList.push_back(makeModExp);
  bkgPdfList.push_back(makeSinglePow);
  bkgPdfList.push_back(makeDoublePow);



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
    bres->SetName("fitres_bkgPdf_"+tag);
    ws->import(*bres);
  }

  ws->import(data);

  return ws;
}


TGraphAsymmErrors* fitChoiceSystematic(RooWorkspace* w, const int nPdfs, const TString * pdfNames, const int *nPars, TString nominal) {
  w->Print();
  RooRealVar * mgg = w->var("mgg");
  std::map<TString, RooAbsPdf*> pdfMap;

  std::map<TString,float> aicTable;
  std::map<TString,float> aicWeightTable;

  //compute the AIC for each pdf and find the minimum
  TString minAICpdf="";
  float minAIC = FLT_MAX;
  for(int i=0;i<nPdfs;i++) {
    std::cout << "fitres_bkgPdf_"+pdfNames[i] << std::endl;
    RooFitResult* res = (RooFitResult*)w->obj("fitres_bkgPdf_"+pdfNames[i]);
    assert(res!=0);
    aicTable[pdfNames[i]] = res->minNll()+2*nPars[i];
    
    if(aicTable[pdfNames[i]] < minAIC) {
      minAIC = aicTable[pdfNames[i]];
      minAICpdf=pdfNames[i];
    }

    //also store the pdfs in the map for ease later
    pdfMap[pdfNames[i]] = w->pdf("bkgPdf_"+pdfNames[i]+"_ext");
    assert(pdfMap[pdfNames[i]]!=0);
  }
  
  //compute the AIC weight for each pdf
  for(int i=0;i<nPdfs;i++) {
    aicWeightTable[pdfNames[i]] = TMath::Exp(-0.5* (aicTable[pdfNames[i]] - minAIC));
  }

  std::map<TString,float> normMap;

  for(int iPdf=0; iPdf<nPdfs; iPdf++) {
    normMap[pdfNames[iPdf]] = pdfMap[pdfNames[iPdf]]->createIntegral(*mgg,RooFit::NormSet(*mgg))->getVal();
  }

  
  float step_size=1.0;

  float mgg_val = mgg->getMin();

  float mgg_max = mgg->getMax();
  int nSteps = (mgg_max-mgg_val)/step_size;

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nSteps);

  int iPt=-1;
  for(; mgg_val <=mgg_max; mgg_val+=step_size) {
    mgg->setVal(mgg_val);
    float nom_v = pdfMap[nominal]->getVal();

    //compute the weighted moments
    float s0=0,s1=0,s2=0;
    for(int iPdf=0; iPdf<nPdfs; iPdf++) {
      float wt = aicWeightTable[pdfNames[iPdf]];
      float v = pdfMap[pdfNames[iPdf]]->getVal()/normMap[pdfNames[iPdf]];
      //std::cout << "\t" << v << "   " << wt << std::endl;
      s0+=wt;
      s1+=wt*v;
      s2+=wt*v*v;
    }

    float mean = s1/s0;
    float std  = TMath::Sqrt(s0*s2-s1*s1)/s0;

    float up = mean+std;
    float down = mean-std;
    //std::cout << mgg_val << " " << nom_v << "   " << mean << " " << up << " " << down << std::endl;

    graph->SetPoint(++iPt,mgg_val,nom_v);
    graph->SetPointError(iPt,0,0,std,std);
  }

  return graph;
}


void computeSigInt(RooWorkspace* w, TString pdfName, TString sfName, float sigEff) {
  RooAbsPdf * pdf = w->pdf("bkgPdf_"+pdfName+"_ext");
  RooRealVar* mgg = w->var("mgg");
  RooFitResult* res = (RooFitResult*)w->obj("fitres_bkgPdf_"+pdfName);

  mgg->setRange("sig",125-2*sigEff, 126+2*sigEff);

  RooAbsReal* integral = pdf->createIntegral(*mgg,RooFit::NormSet(*mgg),RooFit::Range("sig"));
  std::cout << integral->getVal() << "  " << integral->getPropagatedError(*res) << std::endl;
  RooRealVar intRV(sfName,"",integral->getVal());
  intRV.setError(integral->getPropagatedError(*res));
  

  w->import(intRV);
}
