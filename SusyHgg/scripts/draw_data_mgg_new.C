#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TObject.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TH1D.h"

#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooArgList.h"


#include <vector>
#include <algorithm>

struct plot_t {
  TGraph data;
  TGraph pdf;
  TGraphAsymmErrors oneSigma;
  TGraphAsymmErrors twoSigma;

  plot_t(int nData,int nPdf):data(nData),pdf(nPdf),oneSigma(nPdf),twoSigma(nPdf) {}
};


plot_t* buildPlot(const RooWorkspace* w,const float min,const float max, const int nSamples = 2156) {
  if(max<=min) return 0;

  RooAbsPdf* nominalPdf = w->pdf("pdf");
  RooFitResult* fitres = (RooFitResult*)w->obj("fitresult_pdf_data");

  const int nStepsData = int(max-min);
  const int nStepsPdf  = int(max-min)*100;

  const float stepSizeData = (max-min)/nStepsData;
  const float stepSizePdf = (max-min)/nStepsPdf;

  plot_t *p = new plot_t(nStepsData,nStepsPdf);

  RooRealVar* mgg = w->var("mgg");
  mgg->setRange("all",min,max);
  mgg->setBins(nStepsData,"all");

  vector<vector<float> > samples(nStepsPdf, vector<float>(nSamples,0));

  vector<float> sum(nStepsPdf,0);
  vector<float> sumsq(nStepsPdf,0);

  TF1 func("dexp","[0]*([3]*exp([1]*x)+(1-[3])*exp([2]*x))");
  //randomize the parameters nSamples times and build a new PDF
  for(int iSample=0; iSample<nSamples; iSample++) {
    RooArgList l = fitres->randomizePars();
    for(int iPar=0;iPar<4;iPar++) {
      func.SetParameter(iPar,((RooRealVar*)l.at(iPar))->getVal());
    }

    double norm = func.Integral(min,max);

    for(int iStep=0; iStep<nStepsPdf; iStep++) {
      float e = func.Eval(min+iStep*stepSizePdf)*stepSizePdf/norm;
      samples.at(iStep).at(iSample)=e;
      sum.at(iStep)+=e;
      sumsq.at(iStep)+=e*e;
    }
  }

  const int oneSigmaEntries = int(nSamples*(1-0.683)/2);
  const int twoSigmaEntries = int(nSamples*(1-0.955)/2);


  //build the pdfs and evaluate the confidence interval at each point
  for(int iStep=0; iStep<nStepsPdf; iStep++) {
    double m = min+iStep*stepSizePdf;
    mgg->setVal(m);
    
    double pdfVal = nominalPdf->getVal()*stepSizePdf;
    
    double mean = sum.at(iStep)/nSamples;
    double std = sqrt(sumsq.at(iStep)/nSamples - mean*mean);


    //sort the sample predictions at each point
    std::sort(samples.at(iStep).begin(),samples.at(iStep).end());
    double med  = samples.at(iStep).at(nSamples/2);

    p->pdf.SetPoint(iStep,m,pdfVal);
    p->oneSigma.SetPoint(iStep,m,pdfVal);
    p->twoSigma.SetPoint(iStep,m,pdfVal);
    
    //p->oneSigma.SetPointError(iStep,0,0,0.473*std/sqrt(nSamples),0.473*std/sqrt(nSamples));
    //p->twoSigma.SetPointError(iStep,0,0,1.695*std/sqrt(nSamples),1.695*std/sqrt(nSamples));
    if(iStep%100==0) std::cout << m << "   " << pdfVal << "  " << med << "  " << samples.at(iStep).at(oneSigmaEntries) << "  " << samples.at(iStep).at(nSamples-oneSigmaEntries) << std::endl;
    p->oneSigma.SetPointError(iStep,0,0,med- samples.at(iStep).at(oneSigmaEntries), med+ samples.at(iStep).at(nSamples-oneSigmaEntries));
    p->twoSigma.SetPointError(iStep,0,0,med- samples.at(iStep).at(twoSigmaEntries), med+ samples.at(iStep).at(nSamples-twoSigmaEntries));
  }
  
  //build the data histogram
  //RooAbsData* d = w->data("data");
  //RooDataHist rdh("tmp","",mgg,*d);
  



  return p;
}




void draw_data_mgg(TString folderName,bool blind=true,float min=103,float max=160) {
  TFile inputFile(folderName+"/data.root");

  const int nCat = 5;
  TString cats[5] = {"HighPt","Hbb","Zbb","HighRes","LowRes"};

  TCanvas cv;

  for(int iCat=0; iCat < nCat; iCat++) {

    RooWorkspace *ws = (RooWorkspace*)inputFile.Get(cats[iCat]+"_mgg_workspace");

    
    plot_t* plot = buildPlot(ws,min,max);

    TCanvas cv;
    plot->pdf.SetLineColor(kBlue);
    plot->pdf.Draw("AL");

    plot->oneSigma.SetFillColor(kYellow);
    plot->twoSigma.SetFillColor(kGreen);

    plot->oneSigma.SetLineColor(kRed);

    plot->twoSigma.Draw("A4");
    plot->oneSigma.Draw("L4");
    plot->pdf.Draw("L");
    
    TString tag = (blind ? "_BLIND" : "");
    cv.SaveAs(folderName+"/figs/TESTmgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".png");
    cv.SaveAs(folderName+"/figs/TESTmgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".pdf");
    cv.SaveAs(folderName+"/figs/TESTmgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".C");
      
  }
  
}
