#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TObject.h"
#include "TNamed.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

#include "TBox.h"

#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooCurve.h"


void draw_data_mgg(TString folderName,bool blind=true,float min=103,float max=160)
{
  TFile inputFile(folderName+"/data.root");
  
  const int nCat = 5;
  TString cats[5] = {"HighPt","Hbb","Zbb","HighRes","LowRes"};

  TCanvas cv;

  for(int iCat=0; iCat < nCat; iCat++) {

    RooWorkspace *ws  = (RooWorkspace*)inputFile.Get(cats[iCat]+"_mgg_workspace");
    RooFitResult* res = (RooFitResult*)ws->obj("fitresult_pdf_data");

    RooRealVar * mass = ws->var("mgg");
    mass->setRange("all",min,max);
    mass->setRange("blind",121,130);
    mass->setRange("low",106,121);
    mass->setRange("high",130,160);

    mass->setUnit("GeV");
    mass->SetTitle("m_{#gamma#gamma}");
    
    RooAbsPdf * pdf = ws->pdf("pdf");
    RooPlot *plot = mass->frame(min,max,max-min);
    plot->SetTitle("");
    
    RooAbsData* data = ws->data("data")->reduce(Form("mgg > %f && mgg < %f",min,max));
    double nTot = data->sumEntries();
    if(blind) data = data->reduce("mgg < 121 || mgg>130");
    double nBlind = data->sumEntries();
    double norm = nTot/nBlind; //normalization for the plot
    
    data->plotOn(plot,RooFit::Invisible());
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::Range("Full"),RooFit::LineWidth(0.1) );
    plot->Print();
    TH1F* h_tmp = new TH1F("h", "_h", 57, 103, 160);
    TH1F* h_tmp2 = new TH1F("h", "_h", 57, 103, 160);
    h_tmp->SetBinErrorOption(TH1::kPoisson);
    h_tmp2->SetBinErrorOption(TH1::kPoisson);
    data->fillHistogram( h_tmp, *mass);
    h_tmp2->SetLineColor(kBlack);
    h_tmp2->SetMarkerColor(kBlack);
    h_tmp2->SetMarkerStyle(20);
    //h_tmp->SetLineWidth(2);
    h_tmp2->SetStats(0);
    h_tmp->SetLineColor(kBlack);
    h_tmp->SetMarkerColor(kBlack);
    h_tmp->SetMarkerStyle(20);
    //h_tmp->SetLineWidth(2);
    h_tmp->SetStats(0);
    h_tmp->SetTitle("");
    for ( int i = 1; i < h_tmp->GetNbinsX(); i++ )
      {
	h_tmp2->SetBinContent(i,h_tmp->GetBinContent(i));
	if( h_tmp->GetBinContent(i) == 0.0 )
	  {
	    h_tmp->SetBinError(i,1.84); 
	    h_tmp->SetBinContent(i,0.000001);
	  }
	else
	  {
	    h_tmp->SetBinError(i,0); 
	  }
	
      }
    h_tmp->Draw();
    cv.SaveAs("test_test.pdf");
    //add the fix error band
    RooCurve* c = plot->getCurve("pdf_Norm[mgg]_Range[Full]_NormRange[Full]");
    const int Nc = c->GetN();
    //TGraphErrors errfix(Nc);
    //TGraphErrors errfix2(Nc);
    TGraphAsymmErrors errfix(Nc);
    TGraphAsymmErrors errfix2(Nc);
    Double_t *x = c->GetX();
    Double_t *y = c->GetY();
    double NtotalFit = ws->var("Nbkg1")->getVal()*ws->var("Nbkg1")->getVal() + ws->var("Nbkg2")->getVal()*ws->var("Nbkg2")->getVal();
    for( int i = 0; i < Nc; i++ )
      {
	errfix.SetPoint(i,x[i],y[i]);
	errfix2.SetPoint(i,x[i],y[i]);
	mass->setVal(x[i]);      
	double shapeErr = pdf->getPropagatedError(*res)*NtotalFit;
	//double totalErr = TMath::Sqrt( shapeErr*shapeErr + y[i] );
	//total normalization error
	double totalErr = TMath::Sqrt( shapeErr*shapeErr + y[i]*y[i]/NtotalFit ); 
	if ( y[i] - totalErr > .0 )
	  {
	    errfix.SetPointError(i, 0, 0, totalErr, totalErr );
	  }
	else
	  {
	    errfix.SetPointError(i, 0, 0, y[i] - 0.01, totalErr );
	  }
	//2sigma
	if ( y[i] -  2.*totalErr > .0 )
	  {
	    errfix2.SetPointError(i, 0, 0, 2.*totalErr,  2.*totalErr );
	  }
	else
	  {
	    errfix2.SetPointError(i, 0, 0, y[i] - 0.01,  2.*totalErr );
	  }
	/*
	std::cout << x[i] << " " << y[i] << " "
		  << " ,pdf get Val: " << pdf->getVal()
		  << " ,pdf get Prop Err: " << pdf->getPropagatedError(*res)*NtotalFit
		  << " stat uncertainty: " << TMath::Sqrt(y[i]) << " Ntot: " << NtotalFit <<  std::endl;
	*/
      }
    errfix.SetFillColor(kYellow);
    errfix2.SetFillColor(kGreen);


    //pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kGreen),RooFit::Range("Full"), RooFit::VisualizeError(*res,2.0,kFALSE));
    //pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kYellow),RooFit::Range("Full"), RooFit::VisualizeError(*res,1.0,kFALSE));
    //pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kGreen),RooFit::Range("Full"), RooFit::VisualizeError(*res,2.0,kTRUE));
    //pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kYellow),RooFit::Range("Full"), RooFit::VisualizeError(*res,1.0,kTRUE));
    plot->addObject(&errfix,"4");
    plot->addObject(&errfix2,"4");
    plot->addObject(&errfix,"4");
    //data->plotOn(plot);
    TBox blindBox(121,plot->GetMinimum()-(plot->GetMaximum()-plot->GetMinimum())*0.015,130,plot->GetMaximum());
    blindBox.SetFillColor(kGray);
    if(blind) {
      plot->addObject(&blindBox);
      pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kGreen),RooFit::Range("Full"), RooFit::VisualizeError(*res,2.0,kTRUE));
      pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::FillColor(kYellow),RooFit::Range("Full"), RooFit::VisualizeError(*res,1.0,kTRUE));
    }
    //plot->addObject(&errfix,"4");
    //data->plotOn(plot);

    //pdf->plotOn(plot,RooFit::Normalization( norm ) );
    //pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::Range("Full"),RooFit::LineWidth(1.5) );
    pdf->plotOn(plot,RooFit::NormRange( "low,high" ),RooFit::Range("Full"), RooFit::LineWidth(1));
    //data->plotOn(plot);
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
    //h_tmp->SetBinErrorOption(TH1::kPoisson);
    h_tmp2->Draw("sameE0");
    h_tmp2->Draw("sameE1");
    h_tmp->Draw("sameE1");
    TString tag = (blind ? "_BLIND" : "");
    cv.SaveAs(folderName+"/figs/mgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".png");
    cv.SaveAs(folderName+"/figs/mgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".pdf");
    cv.SaveAs(folderName+"/figs/mgg_data_"+cats[iCat]+tag+TString(Form("_%0.0f_%0.0f",min,max))+".C");
      
  }
  
}
