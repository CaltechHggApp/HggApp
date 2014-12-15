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

#include <iostream>
#include <cmath>

#include "getTheoXSec.C"

void makeGraph(TGraphAsymmErrors& graph, TString combine_dir,TString type,bool S0=false);

void compare1DLimit(TString combine_dir_1,TString combine_dir_2,TString lbl1,TString lbl2, TString outputFileName, TString type= "WH",bool S0=true){
  //TString combine_dir = "test_runSusyHgg/signalInj_sms_ChiWH_0_175/";

  TGraphAsymmErrors graph1( (200-125)/25 );
  TGraphAsymmErrors error( (200-125)/25 );
  TGraphAsymmErrors graph2( (200-125)/25 );
  
  TGraphErrors* theo = 0;
  if(type=="WH") theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/CharginoNeutralino.txt");
  //else theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/Higgsino_ElectroHiggs.txt");
  else theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/Higgsino.txt");

  makeGraph(graph1,combine_dir_1,type,S0);
  makeGraph(graph2,combine_dir_2,type,S0);

  TCanvas cv;
  cv.SetLogy(0);

  if(S0) {
    double x,y;
    for(int i=0;i<graph1.GetN();i++) {
      graph1.GetPoint(i,x,y);
      error.SetPoint(i,x,y);
      error.SetPointError(i,0,0, 
			   sqrt( pow(graph1.GetErrorYlow(i),2) + pow(graph2.GetErrorYlow(i),2) ),
			   sqrt( pow(graph1.GetErrorYhigh(i),2) + pow(graph2.GetErrorYhigh(i),2) ) );
			   
      graph1.SetPointError(i,0,0,0,0);
      graph2.SetPointError(i,0,0,0,0);
    }
  }

  graph1.SetMaximum(4);
  graph1.SetMinimum(1e0);
  graph1.GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
  graph1.GetXaxis()->SetTitle("m_{chargino}");
  graph1.SetLineWidth(0.09);
  graph1.SetLineColor(kBlack);
  graph1.SetFillColor(kGreen);
  graph1.Draw("AC");

  if(S0) {
    error.SetFillColor(kGreen);
    error.Draw("3");
    graph1.Draw("C");
  }

  theo->SetFillColor(kBlue);
  theo->SetLineStyle(kDotted);
  theo->SetLineWidth(1);
  theo->Draw("3C");

  graph1.SetLineWidth(0.09);
  graph1.SetLineColor(kBlack);
  graph1.Draw("C");

  graph2.SetLineWidth(0.09);
  graph2.SetLineColor(kRed);
  graph2.Draw("C");

  TLegend leg(0.6,0.7,0.85,0.85);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  if(S0) leg.SetHeader("Expected Limit -- No Systematics");
  else leg.SetHeader("Expected Limit");
  leg.AddEntry(&graph1,lbl1,"l");
  leg.AddEntry(&graph2,lbl2,"l");
  if(S0) leg.AddEntry(&error,"statistical error","f");
  leg.AddEntry(theo,"theoretical","f");

  leg.Draw("SAME");
  
  TLatex prelim(0.65,0.96,"CMS Preliminary");
  prelim.SetNDC();
  prelim.SetTextSize(0.045);
  prelim.Draw();
  
  TLatex lbl(0.5,0.86,"#sqrt{s} = 8 TeV  #int L dt = 19.78 fb^{-1}");
  lbl.SetNDC();
  lbl.SetTextSize(0.045);
  lbl.Draw();
  
  cv.SaveAs(outputFileName);
}


void makeGraph(TGraphAsymmErrors& graph,TString combine_dir,TString type,bool S0) {
  for(int m=125;m<=200;m+=25) {
    int i=(m-125)/25;
    TString file_name = "";
    if(S0) file_name = Form("%s/higgsCombineChi%s_S0_0_%d.Asymptotic.mH120.root",combine_dir.Data(),type.Data(),m);
    else   file_name = Form("%s/higgsCombineChi%s_0_%d.Asymptotic.mH120.root",combine_dir.Data(),type.Data(),m);
    TFile limit_file(file_name);
    TTree *limit_tree = (TTree*)limit_file.Get("limit");
    TTreeFormula limit_form("get_limit","limit",limit_tree);
    
    limit_tree->GetEntry(1);
    float down = limit_form.EvalInstance();
    limit_tree->GetEntry(2); 
    float exp = limit_form.EvalInstance();
    limit_tree->GetEntry(3);
    float up = limit_form.EvalInstance();

    if(i==0) m+=5; //first point is actually at m=130
    graph.SetPoint(i,float(m),exp);
    graph.SetPointError(i,0,0,exp-down,up-exp);
    std::cout << m << ":  " << exp << std::endl;
    if(i==0) m-=5;

    limit_file.Close();
  }
}
