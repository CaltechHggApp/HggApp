#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TObject.h"
#include "TNamed.h"

#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"


void draw_data_mgg(TString folderName,bool blind=true,float min=103,float max=160) {
  TFile inputFile(folderName+"/data.root");

  const int nCat = 5;
  TString cats[5] = {"HighPt","Hbb","Zbb","HighRes","LowRes"};

  TCanvas cv;

  for(int iCat=0; iCat < nCat; iCat++) {

    RooWorkspace *ws = (RooWorkspace*)inputFile.Get(cats[iCat]+"_mgg_workspace");
    RooFitResult* res = (RooFitResult*)ws->obj("fitresult_pdf_data");

    RooRealVar * mass = ws->var("mgg");
    mass->setRange("all",min,max);
    mass->setRange("blind",121,130);
    mass->setRange("low",106,121);
    mass->setRange("high",130,160);

    mass->setUnit("GeV");

    RooAbsPdf * pdf = ws->pdf("pdf");

    double norm=1; //normalization for the plot
    if(blind) { //need a different normalization if we are blinding

      double total = pdf->createIntegral(*mass,RooFit::NormSet(*mass))->getVal();
      double sig   = pdf->createIntegral(*mass,RooFit::NormSet(*mass),RooFit::Range("blind"))->getVal();

      norm = (total-sig)/total;
    }

    mass->SetTitle("m_{#gamma#gamma}");
    RooPlot *plot = mass->frame(min,max,max-min);
    plot->SetTitle("");

    RooAbsData* data = ws->data("data");
    if(blind) data = data->reduce("mgg < 121 || mgg>130");

    data->plotOn(plot);

    pdf->plotOn(plot,RooFit::NormRange( "low,high" ) );
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kGreen), RooFit::VisualizeError(*res,2.0));
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kYellow), RooFit::VisualizeError(*res,1.0));
    data->plotOn(plot);
    //pdf->plotOn(plot,RooFit::Normalization( norm ) );
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ) );

    TLatex lbl0(0.1,0.96,"CMS Preliminary");
    lbl0.SetNDC();
    lbl0.SetTextSize(0.042);
    plot->addObject(&lbl0);
    
    TLatex lbl(0.4,0.96,Form("%s Box",cats[iCat].Data()));
    lbl.SetNDC();
    lbl.SetTextSize(0.042);
    plot->addObject(&lbl);

    TLatex lbl2(0.6,0.96,"#sqrt{s}=8 TeV  L = 19.78 fb^{-1}");
    lbl2.SetNDC();
    lbl2.SetTextSize(0.042);
    plot->addObject(&lbl2);


    int iObj=-1;
    TNamed *obj;
    while( (obj = (TNamed*)plot->getObject(++iObj)) ) {
      obj->SetName(Form("Object_%d",iObj));
    }

    plot->Draw();
    TString tag = (blind ? "_BLIND" : "");
    cv.SaveAs(folderName+"/figs/mgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".png");
    cv.SaveAs(folderName+"/figs/mgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".pdf");
    cv.SaveAs(folderName+"/figs/mgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".C");
      
  }
  
}
