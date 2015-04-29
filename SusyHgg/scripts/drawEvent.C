#include "TLorentzVector.h"
#include "TLine.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TLatex.h"

#include <vector>

void decorateLine(TLine* l) {
  l->SetLineWidth(2);
}

TCanvas * drawEventPhi(TLorentzVector Higgs, std::vector<TLorentzVector>& jets, std::vector<float>& csv, TLorentzVector h1, TLorentzVector h2,TLorentzVector MET,float maxPt=400) {

  TCanvas *c = new TCanvas("event_phi","");
  c->SetGridx(false);
  c->SetGridy(false);

  TH1D *frame_p = new TH1D("frame_p","",2,-1*maxPt,maxPt*1.5);
  frame_p->SetMinimum(-1*maxPt);
  frame_p->SetMaximum(maxPt);

  frame_p->SetXTitle("p_{X} (GeV)");
  frame_p->SetYTitle("p_{Y} (GeV)");

  frame_p->Draw();

  //draw the pt circles

  for(float pt=50; pt<=maxPt; pt+=50) {
    TEllipse* ptCirc = new TEllipse(0,0,pt);
    ptCirc->SetLineColor(1);
    ptCirc->SetLineStyle(3);
    ptCirc->SetFillStyle(0);
    ptCirc->Draw();
  }

  TLine * hLine_p = new TLine(0,0,Higgs.Pt()*TMath::Sin(Higgs.Phi()),Higgs.Pt()*TMath::Cos(Higgs.Phi()));
  hLine_p->SetLineColor(kRed);
  decorateLine(hLine_p);

  hLine_p->Draw();

  TLine* hold_jLine_p=0;
  for(std::vector<TLorentzVector>::const_iterator iJ= jets.begin(); iJ!=jets.end(); iJ++) {
    TLine *jLine_p = new TLine(0,0,iJ->Pt()*TMath::Sin(iJ->Phi()),iJ->Pt()*TMath::Cos(iJ->Phi()));
    jLine_p->SetLineColor(kBlack);
    decorateLine(jLine_p);
    //jLine_p->SetLineStyle(1+iJ-jets.begin());
    jLine_p->Draw();
    if(hold_jLine_p==0) hold_jLine_p=jLine_p;
    if(csv[iJ-jets.begin()]>0.244) {
      TLatex *l = new TLatex(iJ->Pt()*TMath::Sin(iJ->Phi()),iJ->Pt()*TMath::Cos(iJ->Phi()),Form("CSV=%0.2f",csv[iJ-jets.begin()]));
      l->SetTextSize(0.030);
      l->Draw();
    }
  }

  TLine * h1Line_p = new TLine(0,0,h1.Pt()*TMath::Sin(h1.Phi()),h1.Pt()*TMath::Cos(h1.Phi()));
  h1Line_p->SetLineColor(kBlue);
  decorateLine(h1Line_p);
  h1Line_p->SetLineStyle(2);

  h1Line_p->Draw();

  TLine * h2Line_p = new TLine(0,0,h2.Pt()*TMath::Sin(h2.Phi()),h2.Pt()*TMath::Cos(h2.Phi()));
  h2Line_p->SetLineColor(kBlue);
  decorateLine(h2Line_p);
  h2Line_p->SetLineStyle(3);

  h2Line_p->Draw();

  TLine * METLine_p = new TLine(0,0,MET.Pt()*TMath::Sin(MET.Phi()),MET.Pt()*TMath::Cos(MET.Phi()));
  METLine_p->SetLineColor(kGreen);
  decorateLine(METLine_p);

  METLine_p->Draw();


  TLegend *leg = new TLegend(0.7,0.5,0.9,0.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hLine_p,"Higgs","l");
  leg->AddEntry(hold_jLine_p,"Jets","l");
  leg->AddEntry(h1Line_p,"Hemispheres","l");
  leg->AddEntry(METLine_p,"MET","l");

  leg->Draw();

  return c;
}

TCanvas * drawEventEta(TLorentzVector Higgs, std::vector<TLorentzVector>& jets, std::vector<float>& csv, TLorentzVector h1, TLorentzVector h2,TLorentzVector MET,float maxP=400) {

  TCanvas *c = new TCanvas("event_phi","");
  c->SetGridx(false);
  c->SetGridy(false);

  TH1D *frame_p = new TH1D("frame_p","",2,-1*maxP,maxP*1.5);
  frame_p->SetMinimum(-1*maxP);
  frame_p->SetMaximum(maxP);

  frame_p->SetXTitle("p_{#parallel} (GeV)");
  frame_p->SetYTitle("p_{#perp} (GeV)");

  frame_p->Draw();

  //draw the pt circles

  for(float p=50; p<=maxP; p+=50) {
    TEllipse* pCirc = new TEllipse(0,0,p);
    pCirc->SetLineColor(1);
    pCirc->SetLineStyle(3);
    pCirc->SetFillStyle(0);
    pCirc->Draw();
  }

  TLine * hLine_p = new TLine(0,0,Higgs.P()*TMath::Cos(Higgs.Theta()),(Higgs.Phi()<0 ? -1:1)*Higgs.P()*TMath::Sin(Higgs.Theta()));
  hLine_p->SetLineColor(kRed);
  decorateLine(hLine_p);

  hLine_p->Draw();

  TLine* hold_jLine_p=0;
  for(std::vector<TLorentzVector>::const_iterator iJ= jets.begin(); iJ!=jets.end(); iJ++) {
    TLine *jLine_p = new TLine(0,0,iJ->P()*TMath::Cos(iJ->Theta()),(iJ->Phi()<0 ? -1:1)*iJ->P()*TMath::Sin(iJ->Theta()));
    jLine_p->SetLineColor(kBlack);
    decorateLine(jLine_p);
    //jLine_p->SetLineStyle(1+iJ-jets.begin());
    jLine_p->Draw();
    if(hold_jLine_p==0) hold_jLine_p=jLine_p;
    if(csv[iJ-jets.begin()]>0.244) {
      TLatex *l = new TLatex(iJ->P()*TMath::Cos(iJ->Theta()),(iJ->Phi()<0 ? -1:1)*iJ->P()*TMath::Sin(iJ->Theta()),Form("CSV=%0.2f",csv[iJ-jets.begin()]));
      l->SetTextSize(0.030);
      l->Draw();
    }
  }

  TLine * h1Line_p = new TLine(0,0,h1.P()*TMath::Cos(h1.Theta()),(h1.Phi()<0 ? -1:1)*h1.P()*TMath::Sin(h1.Theta()));
  h1Line_p->SetLineColor(kBlue);
  decorateLine(h1Line_p);
  h1Line_p->SetLineStyle(2);
  h1Line_p->Draw();

  TLine * h2Line_p = new TLine(0,0,h2.P()*TMath::Cos(h2.Theta()),(h2.Phi()<0 ? -1:1)*h2.P()*TMath::Sin(h2.Theta()));
  h2Line_p->SetLineColor(kBlue);
  decorateLine(h2Line_p);

  h2Line_p->SetLineStyle(3);
  h2Line_p->Draw();

  TLine * METLine_p = new TLine(0,0,MET.P()*TMath::Cos(MET.Theta()),(MET.Phi()<0 ? -1:1)*MET.P()*TMath::Sin(MET.Theta()));
  METLine_p->SetLineColor(kGreen);
  decorateLine(METLine_p);

  METLine_p->Draw();


  TLegend *leg = new TLegend(0.7,0.5,0.9,0.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hLine_p,"Higgs","l");
  leg->AddEntry(hold_jLine_p,"Jets","l");
  leg->AddEntry(h1Line_p,"Hemispheres","l");
  leg->AddEntry(METLine_p,"MET","l");

  leg->Draw();

  return c;
}


TCanvas * drawEventFromTree(TTree* tree, int iEvt,float maxPt=400,bool eta=false) {
  float pho1_pt,pho1_eta,pho1_phi;
  float pho2_pt,pho2_eta,pho2_phi;

  float ptgg,etagg,phigg,mgg;
  
  float hem1_pt,hem1_eta,hem1_phi,hem1_M;
  float hem2_pt,hem2_eta,hem2_phi,hem2_M;

  int nJets;
  float jet_pt[20],jet_eta[20],jet_phi[20];

  float MET,MET_phi;

  float high_csv,high_csv_pt;
  float sec_csv,sec_csv_pt;

  tree->SetBranchAddress("pho1_pt",&pho1_pt);
  tree->SetBranchAddress("pho1_eta",&pho1_eta);
  tree->SetBranchAddress("pho1_phi",&pho1_phi);

  tree->SetBranchAddress("pho2_pt",&pho2_pt);
  tree->SetBranchAddress("pho2_eta",&pho2_eta);
  tree->SetBranchAddress("pho2_phi",&pho2_phi);

  tree->SetBranchAddress("ptgg",&ptgg);
  tree->SetBranchAddress("etagg",&etagg);
  tree->SetBranchAddress("phigg",&phigg);
  tree->SetBranchAddress("mgg",&mgg);

  tree->SetBranchAddress("hem1_pt",&hem1_pt);
  tree->SetBranchAddress("hem1_eta",&hem1_eta);
  tree->SetBranchAddress("hem1_phi",&hem1_phi);
  tree->SetBranchAddress("hem1_M",&hem1_M);

  tree->SetBranchAddress("hem2_pt",&hem2_pt);
  tree->SetBranchAddress("hem2_eta",&hem2_eta);
  tree->SetBranchAddress("hem2_phi",&hem2_phi);
  tree->SetBranchAddress("hem2_M",&hem2_M);

  tree->SetBranchAddress("t1MET",&MET);
  tree->SetBranchAddress("t1METPhi",&MET_phi);

  tree->SetBranchAddress("Njets",&nJets);
  tree->SetBranchAddress("ptJet",jet_pt);
  tree->SetBranchAddress("etaJet",jet_eta);
  tree->SetBranchAddress("phiJet",jet_phi);

  tree->SetBranchAddress("highest_csv",&high_csv);
  tree->SetBranchAddress("highest_csv_pt",&high_csv_pt);
  
  tree->SetBranchAddress("second_csv",&sec_csv);
  tree->SetBranchAddress("second_csv_pt",&sec_csv_pt);
  
  tree->GetEntry(iEvt);

  TLorentzVector Higgs(0,0,0,0);
  Higgs.SetPtEtaPhiM(ptgg,etagg,phigg,mgg);

  TLorentzVector hem1(0,0,0,0);
  hem1.SetPtEtaPhiM(hem1_pt,hem1_eta,hem1_phi,hem1_M);

  TLorentzVector hem2(0,0,0,0);
  hem2.SetPtEtaPhiM(hem2_pt,hem2_eta,hem2_phi,hem2_M);

  TLorentzVector met(0,0,0,0);
  met.SetPtEtaPhiM(MET,0,MET_phi,0);

  std::vector<TLorentzVector> jets;
  std::vector<float> csv;
  for(int i=0; i<nJets; i++) {
    jets.push_back(TLorentzVector(0,0,0,0));
    jets.back().SetPtEtaPhiM(jet_pt[i],jet_eta[i],jet_phi[i],0);
    if( fabs(high_csv_pt-jet_pt[i])<0.5 ) csv.push_back(high_csv);
    else if( fabs(sec_csv_pt-jet_pt[i])<0.5 ) csv.push_back(sec_csv);
    else csv.push_back(0);
  }

  if(eta)  return drawEventEta(Higgs,jets,csv,hem1,hem2,met,maxPt);
  return drawEventPhi(Higgs,jets,csv,hem1,hem2,met,maxPt);
}
