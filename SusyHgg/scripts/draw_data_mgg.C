#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TObject.h"
#include "TNamed.h"
#include "TGraphErrors.h"

#include "TBox.h"

#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooCurve.h"


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


    mass->SetTitle("m_{#gamma#gamma}");
    RooPlot *plot = mass->frame(min,max,max-min);
    plot->SetTitle("");

    RooAbsData* data = ws->data("data")->reduce(Form("mgg > %f && mgg < %f",min,max));
    double nTot= data->sumEntries();
    if(blind) data = data->reduce("mgg < 121 || mgg>130");
    double nBlind= data->sumEntries();

    double norm=nTot/nBlind; //normalization for the plot


    data->plotOn(plot);

    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::Range("Full"),RooFit::LineWidth(0.1) );

    plot->Print();
    //add the fix error band
    RooCurve* c = plot->getCurve("pdf_Norm[mgg]_Range[Full]_NormRange[Full]");
    const int Nc = c->GetN();
    TGraphErrors errfix(Nc);
    Double_t *x = c->GetX();
    Double_t *y = c->GetY();
      
    for(int i=0;i<Nc;i++) {
      errfix.SetPoint(i,x[i],y[i]);
      mass->setVal(x[i]);      
      //      float relErr = pdf->getPropagatedError(*res)/pdf->getVal();
      //if(relErr<10e-3) relErr = 10e-3;
      float relErr = 1e-2;
      errfix.SetPointError(i,0,relErr*y[i]);
      std::cout << x[i] << " " << y[i] << " " << relErr*y[i] << std::endl;
    }
    errfix.SetFillColor(kYellow);
    


    //pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kGreen),RooFit::Range("Full"), RooFit::VisualizeError(*res,2.0,kFALSE));
    //pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kYellow),RooFit::Range("Full"), RooFit::VisualizeError(*res,1.0,kFALSE));
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kGreen),RooFit::Range("Full"), RooFit::VisualizeError(*res,2.0,kTRUE));
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kYellow),RooFit::Range("Full"), RooFit::VisualizeError(*res,1.0,kTRUE));
    //plot->addObject(&errfix,"4");
    //data->plotOn(plot);
    TBox blindBox(121,plot->GetMinimum()-(plot->GetMaximum()-plot->GetMinimum())*0.015,130,plot->GetMaximum());
    blindBox.SetFillColor(kGray);
    if(blind) {
      plot->addObject(&blindBox);
      pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kGreen),RooFit::Range("Full"), RooFit::VisualizeError(*res,2.0,kTRUE));
      pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kYellow),RooFit::Range("Full"), RooFit::VisualizeError(*res,1.0,kTRUE));
    }
    plot->addObject(&errfix,"4");
    data->plotOn(plot);

    //pdf->plotOn(plot,RooFit::Normalization( norm ) );
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::Range("Full"),RooFit::LineWidth(0.1) );
    
    /*
    pdf->plotOn(plot,RooFit::Normalization(norm),RooFit::Range("all"),RooFit::LineWidth(0.8) );
    //pdf->plotOn(plot,RooFit::Normalization(norm),RooFit::FillColor(kGreen),RooFit::Range("all"), RooFit::VisualizeError(*res,2.0,kFALSE));
    //pdf->plotOn(plot,RooFit::Normalization(norm),RooFit::FillColor(kYellow),RooFit::Range("all"), RooFit::VisualizeError(*res,1.0,kFALSE));
    pdf->plotOn(plot,RooFit::Normalization(norm),RooFit::FillColor(kGreen),RooFit::Range("all"), RooFit::VisualizeError(*res,2.0,kTRUE));
    pdf->plotOn(plot,RooFit::Normalization(norm),RooFit::FillColor(kYellow),RooFit::Range("all"), RooFit::VisualizeError(*res,1.0,kTRUE));
    data->plotOn(plot);
    pdf->plotOn(plot,RooFit::Normalization(norm),RooFit::Range("all"),RooFit::LineWidth(0.8) );
    */
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
