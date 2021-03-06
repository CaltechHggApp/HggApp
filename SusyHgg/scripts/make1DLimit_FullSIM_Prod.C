#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TAxis.h"
#include "TLine.h"

#include "getTheoXSec.C"

void make1DLimit(TString combine_dir,TString type= "WH",bool blind=true){
  //TString combine_dir = "test_runSusyHgg/signalInj_sms_ChiWH_0_175/";

  TGraph obser( (200-125)/25 );
  TGraph graph( (200-125)/25 );
  TGraphAsymmErrors error( (200-125)/25 );

  TGraph obser_r( (200-125)/25 );
  TGraph graph_r( (200-125)/25 );
  TGraphAsymmErrors error_r( (200-125)/25 );
  
  TGraphErrors* theo = 0;
  if(type=="WH") theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/CharginoNeutralino.txt");
  else theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/Higgsino_ElectroHiggs.txt");
  //else theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/Higgsino.txt");

  for(int m=125;m<=200;m+=25) {
    int i=(m-125)/25;
    TFile limit_file(Form("%s/higgsCombineChi%s_0_%d.Asymptotic.mH120.root",combine_dir.Data(),type.Data(),m) );
    TTree *limit_tree = (TTree*)limit_file.Get("limit");
    TTreeFormula limit_form("get_limit","limit",limit_tree);
    
    limit_tree->GetEntry(1);
    float down = limit_form.EvalInstance();
    limit_tree->GetEntry(2);
    float exp = limit_form.EvalInstance();
    limit_tree->GetEntry(3);
    float up = limit_form.EvalInstance();
    limit_tree->GetEntry(5);
    float obs = limit_form.EvalInstance();
    if(i==0) m+=5; //first point is actually at m=130
    graph.SetPoint(i,float(m),exp);
    error.SetPoint(i,float(m),exp);
    error.SetPointError(i,0,0,exp-down,up-exp);
  
    graph_r.SetPoint(i,float(m),exp/theo->Eval(m));
    error_r.SetPoint(i,float(m),exp/theo->Eval(m));
    error_r.SetPointError(i,0,0,(exp-down)/theo->Eval(m),(up-exp)/theo->Eval(m));
  
    obser.SetPoint(i,float(m),obs);
    obser_r.SetPoint(i,float(m),obs/theo->Eval(m));
    if(i==0) m-=5;
  }


    TCanvas cv;
    cv.SetLogy();

    theo->SetMaximum(1e2);
    theo->SetMinimum(1e-2);
    theo->GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
    theo->GetXaxis()->SetTitle("m_{chargino}");
    theo->SetFillColor(kBlue);
    theo->SetLineStyle(kDotted);
    theo->SetLineWidth(1);
    
    error.SetMaximum(1e2);
    error.SetMinimum(1e-2);
    error.GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
    error.GetXaxis()->SetTitle("m_{chargino}");
    error.SetFillColor(kGreen);
    error.Draw("A3");

    theo->Draw("3C");

    graph.SetLineStyle(kDashed);
    graph.SetLineWidth(0.09);
    graph.Draw("C");

    obser.SetLineStyle(1);
    obser.SetLineWidth(1.2);
    if(!blind) obser.Draw("C");

    TLegend leg(0.7,0.7,0.85,0.85);
    leg.SetFillColor(0);
    leg.SetBorderSize(0);
    leg.AddEntry(&graph,"expected","l");
    leg.AddEntry(&error,"expected #pm1#sigma","F");
    leg.AddEntry(theo,"theoretical","f");
    if(!blind)     leg.AddEntry(&obser,"observed","l");

    leg.Draw("SAME");

    TLatex prelim(0.65,0.96,"CMS Preliminary");
    prelim.SetNDC();
    prelim.SetTextSize(0.045);
    prelim.Draw();

    TLatex lbl(0.5,0.86,"#sqrt{s} = 8 TeV  #int L dt = 19.78 fb^{-1}");
    lbl.SetNDC();
    lbl.SetTextSize(0.045);
    lbl.Draw();

    TString infix=(blind ? "" : "_OBS");
    cv.SaveAs(combine_dir+"/expected_exclusion_"+type+"_1D"+infix+".png");
    cv.SaveAs(combine_dir+"/expected_exclusion_"+type+"_1D"+infix+".pdf");

    error_r.SetMaximum(1e2);
    error_r.SetMinimum(1e-2);
    error_r.GetYaxis()->SetTitle("#sigma_{95%}/#sigma_{NLO}");
    if(type=="HH") error_r.GetXaxis()->SetTitle("M_{#chi_{2}^{0}} [GeV]");
    else error_r.GetXaxis()->SetTitle("M_{#chi_{1}^{#pm}} = M_{#chi_{2}^{0}} [GeV]");
    error_r.SetFillColor(kGreen);
    error_r.Draw("A3");

    graph_r.SetLineStyle(kDashed);
    graph_r.SetLineWidth(0.09);
    graph_r.Draw("C");

    TLine l(125,1,205,1);
    l.SetLineWidth(0.18);
    l.Draw("SAME");

    if(!blind) obser_r.Draw("C");
    leg.Draw("SAME");

    lbl.SetY(0.20);
    leg.SetY1NDC(0.28);
    leg.SetY2NDC(0.43);

    prelim.Draw();
    lbl.Draw();
    
    cv.SaveAs(combine_dir+"/expected_exclusion_ratio_"+type+"_1D"+infix+".png");
    cv.SaveAs(combine_dir+"/expected_exclusion_ratio_"+type+"_1D"+infix+".pdf");
}


void make1DLimitHH(TString combine_dir,bool blind=true){
  //TString combine_dir = "test_runSusyHgg/signalInj_sms_ChiHH_0_175/";

  TGraph obser( (275-125)/25 );
  TGraph graph( (275-125)/25 );
  TGraphAsymmErrors error( (275-125)/25 );

    for(int m=125;m<501;m+=25) {
      int i=(m-125)/25;
      TFile limit_file(Form("%s/higgsCombineChiHH_0_%d.Asymptotic.mH120.root",combine_dir.Data(),m) );
      TTree *limit_tree = (TTree*)limit_file.Get("limit");
      TTreeFormula limit_form("get_limit","limit",limit_tree);

      limit_tree->GetEntry(1);
      float down = limit_form.EvalInstance();
      limit_tree->GetEntry(2);
      float exp = limit_form.EvalInstance();
      limit_tree->GetEntry(3);
      float up = limit_form.EvalInstance();
      limit_tree->GetEntry(5);
      float obs = limit_form.EvalInstance();

      graph.SetPoint(i,float(m),exp);
      error.SetPoint(i,float(m),exp);
      error.SetPointError(i,0,0,exp-down,up-exp);

      obser.SetPoint(i,float(m),obs);
    }

    TGraphErrors* theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/Higgsino.txt");

    TCanvas cv;
    cv.SetLogy();

    theo->SetMaximum(1e2);
    theo->SetMinimum(1e-2);
    theo->GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
    theo->GetXaxis()->SetTitle("m_{chargino}");
    theo->SetFillColor(kBlue);
    theo->SetLineStyle(kDotted);
    theo->SetLineWidth(1);
    
    error.SetMaximum(1e2);
    error.SetMinimum(1e-2);
    error.GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
    error.GetXaxis()->SetTitle("m_{chargino}");
    error.SetFillColor(kGreen);
    error.Draw("A3");

    theo->Draw("3C");

    graph.SetLineStyle(kDashed);
    graph.SetLineWidth(0.09);
    graph.Draw("C");

    obser.SetLineStyle(1);
    obser.SetLineWidth(1.2);
    if(!blind) obser.Draw("C");

    TLegend leg(0.7,0.7,0.85,0.85);
    leg.SetFillColor(0);
    leg.SetBorderSize(0);
    leg.AddEntry(&graph,"expected","l");
    leg.AddEntry(&error,"expected #pm1#sigma","F");
    leg.AddEntry(theo,"theoretical","f");
    if(!blind)     leg.AddEntry(&obser,"observed","l");

    leg.Draw("SAME");

    TLatex prelim(0.65,0.96,"CMS Preliminary");
    prelim.SetNDC();
    prelim.SetTextSize(0.045);
    prelim.Draw();

    TLatex lbl(0.5,0.86,"#sqrt{s} = 8 TeV  #int L dt = 19.78 fb^{-1}");
    lbl.SetNDC();
    lbl.SetTextSize(0.045);
    lbl.Draw();

    cv.SaveAs(combine_dir+"expected_exclusion_HH_1D.png");


}
