#include "TString.h"
#include "TGraph.h"
#include "TMultiGraph.h"
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
#include <vector>

#include "getTheoXSec.C"

float lumix = 0.905;
float lumiy = 0.915;
float lumifont = 42;

float cmsx = 0.24;
float cmsy = 0.825;
TString CMSText = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold
float extrax = cmsx + 0.055;
float extray = cmsy - 0.04;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics
// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;
float cmsSize = 0.05;
TString lumiText = "19.8 fb^{-1} (8TeV)";

double wh_125[] = {1.0, 2.16572, 2.90806, 3.99325, 5.52304, 3.05675};
double wh_150[] = {1.17227, 1.63591, 2.26164, 3.16216, 4.23118, 2.43028};
double wh_175[] = {1.05631, 1.25211, 1.78386, 2.60485, 3.39547, 2.0757};
double wh_200[] = {0.776639, 0.898966, 1.39576, 1.94684, 2.66083, 1.51741};
std::vector< double* > wh_limits;

double hh_125[] = {2.48292, 3.6002, 4.91072, 7.71458, 9.87721, 3.8566};
double hh_150[] = {2.31804, 3.3114, 4.40756, 6.03011, 7.99697, 3.63314};
double hh_175[] = {1.30316, 1.61937, 2.32037, 3.3056, 4.77699, 1.91186};
double hh_200[] = {1.08715, 1.17971, 1.78538, 2.46758, 3.2614, 1.66414};
std::vector< double* > hh_limits;

void make1DLimit(TString combine_dir,TString type= "WH",bool blind=true){
  //TString combine_dir = "test_runSusyHgg/signalInj_sms_ChiWH_0_175/";
  //WH
  wh_limits.push_back(wh_125);
  wh_limits.push_back(wh_150);
  wh_limits.push_back(wh_175);
  wh_limits.push_back(wh_200);
  //HH
  hh_limits.push_back(hh_125);
  hh_limits.push_back(hh_150);
  hh_limits.push_back(hh_175);
  hh_limits.push_back(hh_200);

  TGraph obser( (200-125)/25 );
  TGraph graph( (200-125)/25 );
  TGraphAsymmErrors error( (200-125)/25 );
  TGraphAsymmErrors error2S( (200-125)/25 );

  TGraph obser_r( (200-125)/25 );
  TGraph graph_r( (200-125)/25 );
  TGraphAsymmErrors error_r( (200-125)/25 );
  TGraphAsymmErrors error_r2S( (200-125)/25 );
  
  TGraphErrors* theo = 0;
  if(type=="WH") theo = getTheoXSec("xsecs/CharginoNeutralino.txt");
  else theo = getTheoXSec("xsecs/Higgsino_ElectroHiggs.txt");
  //else theo = getTheoXSec("/home/amott/HggApp/SusyHgg/xsecs/Higgsino.txt");

  for(int m=125;m<=200;m+=25) {
    int i=(m-125)/25;
    TFile limit_file(Form("%s/higgsCombineChi%s_0_%d.Asymptotic.mH120.root",combine_dir.Data(),type.Data(),m) );
    TTree *limit_tree = (TTree*)limit_file.Get("limit");
    TTreeFormula limit_form("get_limit","limit",limit_tree);
    
    float down_2s = -1;
    float down = -1;
    float exp = -1;
    float up = -1;
    float up_2s = -1;
    float obs = -1;

    if( type == "WH" )
      {
	down_2s = wh_limits.at(i)[0];
	down = wh_limits.at(i)[1];
	exp = wh_limits.at(i)[2];
	up = wh_limits.at(i)[3];
	up_2s = wh_limits.at(i)[4];
	obs = wh_limits.at(i)[5];
      }
    else if ( type == "HH")
      {
	down_2s = hh_limits.at(i)[0];
	down = hh_limits.at(i)[1];
	exp = hh_limits.at(i)[2];
	up = hh_limits.at(i)[3];
	up_2s = hh_limits.at(i)[4];
	obs = hh_limits.at(i)[5];
      }
    else
      {
	std::cerr << "UNRECOGNIZED OPTION!!! QUITTING" << std::endl;
      }
    
    if(i==0) m+=5; //first point is actually at m=130
    graph.SetPoint(i,float(m), exp);
    error.SetPoint(i,float(m), exp);
    error2S.SetPoint(i, float(m), exp);
    error.SetPointError(i, 0, 0, exp-down, up-exp);
    error2S.SetPointError(i, 0 , 0 , exp-down_2s, up_2s-exp);
  
    graph_r.SetPoint(i,float(m),exp/theo->Eval(m));
    error_r.SetPoint(i,float(m),exp/theo->Eval(m));
    error_r2S.SetPoint(i,float(m),exp/theo->Eval(m));
    error_r.SetPointError(i,0,0,(exp-down)/theo->Eval(m),(up-exp)/theo->Eval(m));
    error_r2S.SetPointError(i, 0, 0, (exp-down_2s)/theo->Eval(m), (up_2s-exp)/theo->Eval(m) );
    
  
    obser.SetPoint(i,float(m),obs);
    obser_r.SetPoint(i,float(m),obs/theo->Eval(m));
    if(i==0) m-=5;
  }


    TCanvas cv;
    cv.SetLogy();
    cv.SetGrid(1,1);
    theo->SetMaximum(1e2);
    theo->SetMinimum(1e-2);
    theo->GetYaxis()->SetLabelSize(0.05);
    theo->GetYaxis()->SetTitleSize(0.06);
    theo->GetYaxis()->SetTitleOffset(0.8);
    theo->GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
    theo->GetXaxis()->SetTitle("m_{chargino} (GeV)");
    if(type=="HH") theo->GetXaxis()->SetTitle("m_{neutralino} (GeV)");
    theo->SetFillColor(kBlue);
    theo->SetLineStyle(kDotted);
    theo->SetLineWidth(2);
    
    error.SetMaximum(1e2);
    error.SetMinimum(1e-2);
    error.GetYaxis()->SetLabelSize(0.04);
    error.GetYaxis()->SetTitleSize(0.06);
    error.GetYaxis()->SetTitleOffset(0.8);
    error.GetXaxis()->SetLabelSize(0.04);
    error.GetXaxis()->SetTitleSize(0.05);
    error.GetXaxis()->SetTitleOffset(0.9);
    error.GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
    error.GetXaxis()->SetTitle("m_{chargino} (GeV)");
    if(type=="HH") error.GetXaxis()->SetTitle("m_{neutralino} (GeV)");
    error.SetFillColor(kGreen);
    error2S.SetFillColor(kYellow);
    error2S.SetTitle("");
    error.SetTitle("");
    error.Draw("A3");
    error2S.Draw("3SAME");
    error.Draw("3");
    
    theo->SetTitle("");
    theo->Draw("3C");

    graph.SetLineStyle(kDashed);
    graph.SetLineWidth(2.0);
    graph.SetTitle("");
    graph.Draw("C");

    obser.SetLineStyle(1);
    obser.SetLineWidth(2.0);
    obser.SetTitle("");
    if(!blind) obser.Draw("C");

    TLegend leg(0.65,0.65,0.89,0.89);
    leg.SetFillColor(0);
    leg.SetBorderSize(0);
    leg.AddEntry(&graph,"expected","l");
    leg.AddEntry(&error,"expected #pm1#sigma","F");
    leg.AddEntry(&error2S,"expected #pm2#sigma","F");
    leg.AddEntry(theo,"theoretical","f");
    if(!blind)     leg.AddEntry(&obser,"observed","l");

    leg.Draw("SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    
    float extraTextSize = extraOverCmsTextSize*cmsSize;
    latex.SetTextFont(lumifont);
    latex.SetTextAlign(31); 
    latex.SetTextSize(cmsSize);    
    latex.DrawLatex(lumix, lumiy,lumiText);
    
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(31); 
    latex.SetTextSize(cmsSize);
    latex.DrawLatex(cmsx, cmsy, CMSText);
    
    latex.SetTextFont(extraTextFont);
    latex.SetTextAlign(31); 
    latex.SetTextSize(extraTextSize);
    latex.DrawLatex(extrax, extray, extraText);

    TString infix=(blind ? "" : "_OBS");
    cv.SaveAs(combine_dir+"/expected_exclusion_"+type+"_1D"+infix+"_v2.png");
    cv.SaveAs(combine_dir+"/expected_exclusion_"+type+"_1D"+infix+"_v2.pdf");
    cv.SaveAs(combine_dir+"/expected_exclusion_"+type+"_1D"+infix+"_v2.C");

    error_r.SetMaximum(1e2);
    error_r.SetMinimum(1e-2);

    error_r.SetTitle("");
    error_r.GetYaxis()->SetLabelSize(0.04);
    error_r.GetYaxis()->SetTitleSize(0.06);
    error_r.GetYaxis()->SetTitleOffset(0.8);
    error_r.GetXaxis()->SetLabelSize(0.04);
    error_r.GetXaxis()->SetTitleSize(0.05);
    error_r.GetXaxis()->SetTitleOffset(0.9);
    error_r.GetYaxis()->SetTitle("#sigma_{95%}/#sigma_{NLO}");
    if(type=="HH") error_r.GetXaxis()->SetTitle("m_{neutralino} (GeV)");
    else error_r.GetXaxis()->SetTitle("m_{chargino} (GeV)");
    error_r.SetFillColor(kGreen);
    error_r2S.SetFillColor(kYellow);
    error_r2S.SetTitle("");
    error_r.Draw("A3");
    error_r2S.Draw("3SAME");
    error_r.Draw("3SAME");

    graph_r.SetLineStyle(kDashed);
    graph_r.SetLineWidth(2);
    graph_r.SetTitle("");
    graph_r.Draw("C");

    TLine l(125,1,205,1);
    l.SetLineWidth(3);
    l.SetLineColor(kBlue);
    //l.Draw("SAME");

    obser_r.SetLineWidth(2);
    if(!blind) obser_r.Draw("C");
    leg.Draw("SAME");

    //lbl.SetY(0.20);
    //leg.SetY1NDC(0.28);
    //leg.SetY2NDC(0.43);

    //prelim.Draw();
    //lbl.Draw();
    latex.SetTextFont(lumifont);
    latex.SetTextAlign(31); 
    latex.SetTextSize(cmsSize);    
    latex.DrawLatex(lumix, lumiy,lumiText);
    
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(31); 
    latex.SetTextSize(cmsSize);
    latex.DrawLatex(cmsx, cmsy, CMSText);
    
    latex.SetTextFont(extraTextFont);
    latex.SetTextAlign(31); 
    latex.SetTextSize(extraTextSize);
    latex.DrawLatex(extrax, extray, extraText);
    
    cv.SaveAs(combine_dir+"/expected_exclusion_ratio_"+type+"_1D"+infix+"_v2.png");
    cv.SaveAs(combine_dir+"/expected_exclusion_ratio_"+type+"_1D"+infix+"_v2.pdf");
    cv.SaveAs(combine_dir+"/expected_exclusion_ratio_"+type+"_1D"+infix+"_v2.C");
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
    theo->SetLineWidth(2.0);
    
    error.SetMaximum(1e2);
    error.SetMinimum(1e-2);
    error.GetYaxis()->SetTitle("95% CL #sigma upper limit (pb)");
    error.GetXaxis()->SetTitle("m_{chargino}");
    error.SetFillColor(kGreen);
    error.Draw("A3");

    theo->Draw("3C");

    graph.SetLineStyle(kDashed);
    graph.SetLineWidth(2);
    graph.Draw("C");

    obser.SetLineStyle(1);
    obser.SetLineWidth(2);
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
