#include "Plotter.hpp"

Plotter::Plotter(TString inputFileName, TString outDir, TString outTag) {
  inputFile = new TFile(inputFileName);
  inputWs   = (RooWorkspace*)inputFile->Get("susy_hgg_workspace");

  outputTag = outTag;
  outputDir = outDir;

  assert(inputWs!=0);

  RooCategory *evtcat = (RooCategory*)inputWs->obj("evtcat");
  getCategories(evtcat);

  mgg = inputWs->var("mgg");
  mgg->setRange("plot",110,140);
  mgg->setBins(20,"plot");
  //mgg->setRange("sigplot",115,135);
  //mgg->setBins(40,"sigplot");

  setStyle();
}

Plotter::~Plotter() {
  inputFile->Close();
}

void Plotter::getCategories(RooCategory* roocat) {
  assert(roocat != 0);
  int i=0;
  while( roocat->setIndex(i++) == 0 ) {
    categories.push_back(roocat->getLabel());
  }
}

void Plotter::plotMggFits(TString catName) {
  RooPlot* frame = mgg->frame(110,140,20);
  inputWs->data("data_Combined")->reduce("evtcat==evtcat::"+catName)->plotOn(frame);
  inputWs->pdf(catName+"_bkgModel")->plotOn(frame);

  TCanvas cv;
  frame->Draw();
  cv.SaveAs( outputDir+"/mgg_data_"+catName+"_"+outputTag+".png");
}

void Plotter::plotSignalPeaks(TString catName) {
  RooPlot *frame = mgg->frame(115,135,40);
  inputWs->data("SMTot_Combined")->reduce("evtcat==evtcat::"+catName)->plotOn(frame);

  TLatex *owner = new TLatex(0.12,0.96,"CMS Preliminary Simulation");
  owner->SetNDC();
  owner->SetTextSize(0.045);
  owner->SetTextColor(kBlack);
  frame->addObject(owner);

  TLatex *process = new TLatex(0.55,0.88,"SM H#rightarrow#gamma#gamma");
  process->SetNDC();
  process->SetTextSize(0.045);
  process->SetTextColor(kBlack);
  frame->addObject(process);


  TLatex *catLbl = new TLatex(0.55,0.8,catName+" Category");
  catLbl->SetNDC();
  catLbl->SetTextSize(0.045);
  catLbl->SetTextColor(kBlack);
  frame->addObject(catLbl);

  RooRealVar *se = inputWs->var(catName+"_SMTot_sigEff");
  TLatex *sigeff = new TLatex(0.55,0.72,Form("#sigma_{eff} = %0.2f GeV",se->getVal()));
  sigeff->SetNDC();
  sigeff->SetTextSize(0.045);
  sigeff->SetTextColor(kBlack);
  frame->addObject(sigeff);

  TCanvas cv;
  frame->Draw();
  cv.SaveAs( outputDir+"/mgg_SMTot_"+catName+"_"+outputTag+".png");  
}

void Plotter::plotTH2s() {
  TList *l = inputFile->GetListOfKeys();

  TCanvas cv;

  for( int i=0; i<l->GetEntries();i++) {
    TObject *o = inputFile->Get(l->At(i)->GetName());
    if( strncmp(o->IsA()->GetName(),"TH2F",4)!=0 ) continue;

    if( strstr(o->GetName(),"SigRegions")==NULL ) o->Draw("COLZ");
    else o->Draw("COLZTEXT");

    if( strstr(o->GetName(),"sidebandSub")==NULL ) cv.SetLogz();
    else cv.SetLogz(false);
    cv.SaveAs( outputDir+"/RMR_"+l->At(i)->GetName()+".png");
  }
}

void Plotter::plotFrenchFlag(TString catName) {
  TH2F* sm = (TH2F*)inputFile->Get( Form("SMTot_%s_Signal_SigRegions",catName.Data()) );
  TH2F* bkg = (TH2F*)inputFile->Get( Form("data_%s_Sideband_SigRegions",catName.Data()) );
  TH2F* obs = (TH2F*)inputFile->Get( Form("data_%s_Signal_SigRegions",catName.Data()) );

  TH2F* ff = (TH2F*)obs->Clone( Form("data_%s_FrenchFlag",catName.Data()) );

  RooRealVar* scale = inputWs->var( Form("bkg_scale_factor_%s",catName.Data()) );
  RooRealVar* sig_int = inputWs->var( Form("signal_fit_int_%s",catName.Data()) );

  for(int iX=1;iX<ff->GetNbinsX()+1;iX++) {
    for(int iY=1;iY<ff->GetNbinsY()+1;iY++) {
      float bkg_tot = sm->GetBinContent(iX,iY)+bkg->GetBinContent(iX,iY)*sig_int->getVal()/bkg->Integral();
      float bkg_err = sqrt( sm->GetBinContent(iX,iY) + 
			    pow( bkg->GetBinContent(iX,iY)*scale->getVal(),2)*( 1/bkg->GetBinContent(iX,iY) + pow(sig_int->getError()/sig_int->getVal(),2)) );

      if(bkg->GetBinContent(iX,iY)==0) bkg_err = sqrt(sm->GetBinContent(iX,iY));
      if(bkg_err<1)bkg_err = 1; //HACK!!

      ff->SetBinContent(iX,iY,(ff->GetBinContent(iX,iY)-bkg_tot)/bkg_err);
    }
  }

  TCanvas cv;

  ff->Draw("COLZTEXT");

  cv.SaveAs( outputDir+"/FrenchFlag_"+catName+".png");

  delete ff;
}

void Plotter::Run() {
  plotTH2s();
  for(auto& it: categories) {
    if(!isSMS) {
      plotMggFits(it);
      plotSignalPeaks(it);
      plotFrenchFlag(it);
    }
  }
}


void Plotter::setStyle(){

  vecbosStyle = new TStyle("vecbosStyle","Style for P-TDR");

  // For the canvas:
  vecbosStyle->SetCanvasBorderMode(0);
  vecbosStyle->SetCanvasColor(kWhite);
  vecbosStyle->SetCanvasDefH(600); //Height of canvas
  vecbosStyle->SetCanvasDefW(900); //Width of canvas
  vecbosStyle->SetCanvasDefX(0);   //POsition on screen
  vecbosStyle->SetCanvasDefY(0);

  // For the Pad:
  vecbosStyle->SetPadBorderMode(0);
  // vecbosStyle->SetPadBorderSize(Width_t size = 1);
  vecbosStyle->SetPadColor(kWhite);
  vecbosStyle->SetPadGridX(true);
  vecbosStyle->SetPadGridY(true);
  vecbosStyle->SetGridColor(0);
  vecbosStyle->SetGridStyle(3);
  vecbosStyle->SetGridWidth(1);

  // For the frame:
  vecbosStyle->SetFrameBorderMode(0);
  vecbosStyle->SetFrameBorderSize(1);
  vecbosStyle->SetFrameFillColor(0);
  vecbosStyle->SetFrameFillStyle(0);
  vecbosStyle->SetFrameLineColor(1);
  vecbosStyle->SetFrameLineStyle(1);
  vecbosStyle->SetFrameLineWidth(1);

  // set the paper & margin sizes
  vecbosStyle->SetPaperSize(20,26);
  vecbosStyle->SetPadTopMargin(0.05);
  vecbosStyle->SetPadRightMargin(0.13);
  vecbosStyle->SetPadBottomMargin(0.16);
  vecbosStyle->SetPadLeftMargin(0.12);

  // use large Times-Roman fonts
  vecbosStyle->SetTitleFont(132,"xyz");  // set the all 3 axes title font
  vecbosStyle->SetTitleFont(132," ");    // set the pad title font
  vecbosStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  vecbosStyle->SetTitleSize(0.06," ");   // set the pad title size
  vecbosStyle->SetLabelFont(132,"xyz");
  vecbosStyle->SetLabelSize(0.05,"xyz");
  vecbosStyle->SetLabelColor(1,"xyz");
  vecbosStyle->SetTextFont(132);
  vecbosStyle->SetTextSize(0.08);
  vecbosStyle->SetStatFont(132);

  vecbosStyle->SetTitleOffset(0.9,"Y");
  // use bold lines and markers
  vecbosStyle->SetMarkerStyle(8);
  vecbosStyle->SetMarkerSize(1.2);
  vecbosStyle->SetHistLineWidth(1.85);
  vecbosStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //..Get rid of X error bars
  //vecbosStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  vecbosStyle->SetOptTitle(0);
  vecbosStyle->SetOptStat(0);
  vecbosStyle->SetOptFit(11111111);

  // put tick marks on top and RHS of plots
  vecbosStyle->SetPadTickX(1);
  vecbosStyle->SetPadTickY(1);

  // set a decent palette
  vecbosStyle->SetPalette(1);

  vecbosStyle->cd();

  gROOT->SetStyle("vecbosStyle");
  gROOT->ForceStyle();

}
