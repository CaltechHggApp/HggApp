#include "MakeSpinPlots.h"
#include "MakeSpinWorkspace.h"
#include "subtract.cc"


MakeSpinPlots::MakeSpinPlots(TString inputFileName, TString outTag, TString workspaceName){
  isSetup=false;
  inputFile = new TFile(inputFileName);
  ws = (RooWorkspace*)inputFile->Get(workspaceName);

  MakeSpinFits::getLabels("fitlabels",&mcNames,ws);

  if(mcNames.size()==0) MakeSpinFits::getLabels("labels",&mcNames,ws);
  MakeSpinFits::getLabels("evtcat",&catNames,ws);
  MakeSpinFits::getLabels("CosThetaBins",&cosThetaBins,ws);

  basePath = "figs/";
  outputTag=outTag;

  smName="";

  setStyle();

}

MakeSpinPlots::~MakeSpinPlots(){
  inputFile->Close();
}

void MakeSpinPlots::runAll(){
  if(!isSetup) setupDir();
  std::vector<TString>::const_iterator mcIt = mcNames.begin();
  for(; mcIt != mcNames.end(); mcIt++){
    runAll(*mcIt);
  }
}

void MakeSpinPlots::runAll(TString mcName){
  if(!isSetup) setupDir();
  if(ws->var(Form("Data_%s_FULLFIT_Nsig",mcName.Data()))==0) return;

  std::vector<TString>::const_iterator catIt = catNames.begin();
  for(; catIt != catNames.end(); catIt++){
    runAll(*catIt,mcName);
  }

  DrawSpinSubTotBackground(mcName,false);
  DrawSpinSubTotBackground(mcName,true);
}

void MakeSpinPlots::runAll(TString tag, TString mcName){
  if(!isSetup) setupDir();
  getFitValues(tag,mcName);
  DrawBlindFit(tag,mcName);
  DrawFit(tag,mcName);
  DrawIndFit(tag,mcName);
  PlotSignalFits(tag,mcName);
  for(std::vector<TString>::const_iterator csBinIt = cosThetaBins.begin();
      csBinIt != cosThetaBins.end(); csBinIt++){
    TString tmp = tag+"_"+*csBinIt;
    double sigEff = ws->var(Form("%s_FIT_%s_sigmaEff",mcName.Data(),tmp.Data()))->getVal();
    fitSigEff[tPair(mcName,tmp)]    = dPair(sigEff,0);
    PlotSignalFits(tag,mcName,*csBinIt);    
    DrawBlindFit(tag,mcName,*csBinIt);
    DrawFit(tag,mcName,*csBinIt);
  }
  DrawSpinBackground(tag,mcName,false);
  DrawSpinBackground(tag,mcName,true);
  DrawSpinSubBackground(tag,mcName,false);
  DrawSpinSubBackground(tag,mcName,true);
}

void MakeSpinPlots::getFitValues(TString tag,TString mcName){
  std::cout << tag << std::endl;
  //Data_RSG125_FULLFIT_Nsig
  double sig  = ws->var(Form("Data_%s_FULLFIT_Nsig",mcName.Data()))->getVal();
  double sige = ws->var(Form("Data_%s_FULLFIT_Nsig",mcName.Data()))->getError();
  double fsig = ws->var(Form("Data_%s_FULLFIT_%s_fsig",mcName.Data(),tag.Data()))->getVal();
  std::cout << sig*fsig << " +- " << sige*fsig <<std::endl;
  //double bkg  = ws->var(Form("Data_%s_FULLFIT_Nbkg",mcName.Data()))->getVal();
  //double bkge = ws->var(Form("Data_%s_FULLFIT_Nbkg",mcName.Data()))->getError();
  double bkg = ws->var(Form("Data_%s_FULLFIT_%s_Nbkg",mcName.Data(),tag.Data()))->getVal();
  double bkge = ws->var(Form("Data_%s_FULLFIT_%s_Nbkg",mcName.Data(),tag.Data()))->getError();
  //std::cout << bkg*fbkg << " +- " << bkge*fbkg <<std::endl;
  double mean = ws->var(Form("%s_FIT_%s_mean",mcName.Data(),tag.Data()))->getVal();
  double meane= ws->var(Form("%s_FIT_%s_mean",mcName.Data(),tag.Data()))->getError();

  std::cout << mean << " +- " << meane <<std::endl;  

  double sigEff = ws->var(Form("%s_FIT_%s_sigmaEff",mcName.Data(),tag.Data()))->getVal();
  double sigEffE = 0;
  std::cout << sigEff << " +- " << sigEffE <<std::endl;  
  //Data_BKGFIT_EB_0_bkgModel
  RooAbsPdf * bkgPdf = ws->pdf(Form("Data_BKGFIT_%s_bkgModel",tag.Data()));
  std::cout << bkgPdf <<std::endl;
  RooRealVar range("range","",mean-sigEff,mean+sigEff);
  RooRealVar all("all","",ws->var("mass")->getMin(),ws->var("mass")->getMax());
  double BkgInRange  = bkgPdf->createIntegral(range)->getVal()/bkgPdf->createIntegral(all)->getVal()*bkg;
  double BkgInRangeE = bkgPdf->createIntegral(range)->getVal()/bkgPdf->createIntegral(all)->getVal()*bkge;
  

  tPair lbl(mcName,tag);

  nSignal[lbl]      = dPair(sig*fsig,sige*fsig);
  nBackground[lbl]  = dPair(bkg,bkge);
  fitMean[lbl]      = dPair(mean,meane);
  fitSigEff[lbl]    = dPair(sigEff,sigEffE);
  fitBkg1Sigma[lbl] = dPair(BkgInRange,BkgInRangeE);

}


void MakeSpinPlots::DrawBlindFit(TString tag, TString mcName,TString cosThetaBin){
  TString fitTag="FULLFIT";
  TString cat = "evtcat";
  if(cosThetaBin!=""){
    tag+="_"+cosThetaBin;
    fitTag="FULLCOSTFIT";
    cat = "evtcat_cosT";
  }
  TString dataTag = "_Combined";
  if(cosThetaBin!="") dataTag+="_CosTBin";

  TCanvas *cv = new TCanvas(Form("%s_%s",mcName.Data(),tag.Data()));

  
  RooRealVar* mass = ws->var("mass");
  mass->setBins( (mass->getMax() - mass->getMin())/1.5 ); //enfore 1.5GeV bin width
  RooPlot* frame  = mass->frame();
  double Nb = ws->var(Form("Data_BKGFIT_%s_Nbkg",tag.Data()))->getVal();
  cout << Nb << endl;
  RooDataSet *blind = (RooDataSet*)ws->data("Data"+dataTag)->reduce(TString("((mass<119) || (mass>135.5)) && ")+cat+"=="+cat+"::"+tag);
  blind->plotOn(frame);

  tPair lbl(mcName,tag);
  double nBkg = ws->data("Data"+dataTag)->sumEntries(cat+"=="+cat+"::"+tag);

  ws->pdf( Form("Data_BKGFIT_%s_bkgModel",tag.Data()) )->plotOn(frame,RooFit::Range("all"),RooFit::Normalization(nBkg/blind->sumEntries()),
									       RooFit::LineColor(kRed));

  //TLatex *prelim = new TLatex(250,x->GetXmax()-40.,"CMS Preliminary");
  TLatex *prelim = new TLatex(0.12,0.96,"CMS Preliminary");
  TLatex *lum = new TLatex(0.7,0.96,Form("#sqrt{s}=8 TeV  L = %0.1f fb^{-1}",lumi));
  prelim->SetNDC();
  lum->SetNDC();
  prelim->SetTextSize(0.045);
  prelim->SetTextColor(kBlack);
  lum->SetTextSize(0.045);
  lum->SetTextColor(kBlack);

  TLatex *owner = new TLatex(0.6,0.88,"Alex Mott (Nov. 13, 2012)");
  owner->SetNDC();
  owner->SetTextSize(0.045);
  owner->SetTextColor(kBlack);

  TLatex *Nbkg = new TLatex(0.7,0.8,Form("N_{bkg}= %0.1f #pm %0.1f",nBackground[lbl].first,nBackground[lbl].second));
  Nbkg->SetNDC();
  Nbkg->SetTextSize(0.045);

  TLatex *sig = new TLatex(0.7,0.72,Form("#sigma_{eff} = %0.1f #pm %0.2f",fitSigEff[lbl].first,fitSigEff[lbl].second));
  sig->SetNDC();
  sig->SetTextSize(0.045);

  TLatex *expBkg = new TLatex(0.7,0.64,Form("B @ 125 = %0.1f",fitBkg1Sigma[lbl].first));
  expBkg->SetNDC();
  expBkg->SetTextSize(0.045);


  frame->addObject(prelim);
  frame->addObject(lum);
  //frame->addObject(owner);
  frame->addObject(Nbkg);
  frame->addObject(sig);
  frame->addObject(expBkg);
  frame->Draw();
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s_BLIND.png",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/C/mgg-%s-%s-%s_BLIND.C",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s_BLIND.pdf",outputTag.Data(),mcName.Data(),tag.Data()) );
  delete cv;
}

void MakeSpinPlots::DrawFit(TString tag, TString mcName,TString cosThetaBin){
  TString fitTag="FULLFIT";
  TString cat = "evtcat";
  if(cosThetaBin!=""){
    tag+="_"+cosThetaBin;
    fitTag="FULLCOSTFIT";
    cat = "evtcat_cosT";
 }
  TString dataTag = "_Combined";
  if(cosThetaBin!="") dataTag+="_CosTBin";
  TCanvas *cv = new TCanvas(Form("%s_%s",mcName.Data(),tag.Data()));

  RooRealVar* mass = ws->var("mass");
  RooPlot* frame  = mass->frame();

  mass->setBins( (mass->getMax() - mass->getMin())/1.5 ); //enfore 1.5GeV bin width

  tPair lbl(mcName,tag);
  std::cout << Form("Data_%s_%s_%s_Nsig",mcName.Data(),fitTag.Data(),tag.Data()) <<std::endl;
  double Ns = ((RooFormulaVar*)ws->obj( Form("Data_%s_%s_%s_Nsig",mcName.Data(),fitTag.Data(),tag.Data()) ))->getVal();
  double Nb = ((RooFormulaVar*)ws->obj( Form("Data_%s_%s_%s_Nbkg",mcName.Data(),fitTag.Data(),tag.Data()) ))->getVal();

  double Nblind = ws->data("Data"+dataTag)->reduce("(mass>100 && mass<119) || (mass>135.5 && mass<170)")->sumEntries(cat+"=="+cat+"::"+tag);
  double Ntot   = ws->data("Data"+dataTag)->sumEntries(cat+"=="+cat+"::"+tag);

  RooFitResult* fitres = (RooFitResult*)ws->obj(Form("Data_%s_%s_fitResult",mcName.Data(),fitTag.Data())); 
  std::cout << fitres << std::endl;
  ws->data("Data"+dataTag)->reduce(cat+"=="+cat+"::"+tag)->plotOn(frame,RooFit::LineColor(kWhite),RooFit::MarkerColor(kWhite));
  //Data_Hgg125_FULLFIT_EB_0
  ws->pdf(Form("Data_%s_%s_%s",mcName.Data(),fitTag.Data(),tag.Data()))->plotOn(frame, RooFit::FillColor(kGreen),RooFit::VisualizeError(*fitres,2.0));
  ws->pdf(Form("Data_%s_%s_%s",mcName.Data(),fitTag.Data(),tag.Data()))->plotOn(frame, RooFit::FillColor(kYellow),RooFit::VisualizeError(*fitres,1.0));
  ws->pdf(Form("Data_%s_%s_%s",mcName.Data(),fitTag.Data(),tag.Data()))->plotOn(frame, RooFit::LineColor(kRed));
  std::cout << "1" << std::endl;
  ws->pdf(Form("Data_BKGFIT_%s_bkgModel",tag.Data()))->plotOn(frame, RooFit::Normalization(Nb/(Nb+Ns)),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  std::cout << "2" << std::endl;

  ws->data("Data"+dataTag)->reduce(cat+"=="+cat+"::"+tag)->plotOn(frame);
  frame->Draw();

  //TLatex *prelim = new TLatex(250,x->GetXmax()-40.,"CMS Preliminary");
  TLatex *prelim = new TLatex(0.12,0.96,"CMS Preliminary");
  TLatex *lum = new TLatex(0.7,0.96,Form("#sqrt{s}=8 TeV  L = %0.1f fb^{-1}",lumi));
  prelim->SetNDC();
  lum->SetNDC();
  prelim->SetTextSize(0.045);
  prelim->SetTextColor(kBlack);
  lum->SetTextSize(0.045);
  lum->SetTextColor(kBlack);

  TLatex *owner = new TLatex(0.6,0.88,"Alex Mott (Nov. 13, 2012)");
  owner->SetNDC();
  owner->SetTextSize(0.045);
  owner->SetTextColor(kBlack);

  TLatex *mu = new TLatex(0.7,0.8,Form("#mu = %0.1f #pm %0.2f", fitMean[lbl].first,fitMean[lbl].second));
  mu->SetNDC();
  mu->SetTextSize(0.045);

  TLatex *sig = new TLatex(0.7,0.72,Form("#sigma_{eff} = %0.1f #pm %0.2f", fitSigEff[lbl].first,fitSigEff[lbl].second));
  sig->SetNDC();
  sig->SetTextSize(0.045);

  TLatex *Nsig = new TLatex(0.7,0.64,Form("N_{sig}= %0.1f #pm %0.1f",nSignal[lbl].first,nSignal[lbl].second));
  Nsig->SetNDC();
  Nsig->SetTextSize(0.045);


  frame->addObject(prelim);
  frame->addObject(lum);
  //frame->addObject(owner);
  frame->addObject(mu);
  frame->addObject(sig);
  frame->addObject(Nsig);
  frame->Draw();
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s.png",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/C/mgg-%s-%s-%s.C",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s.pdf",outputTag.Data(),mcName.Data(),tag.Data()) );
  delete cv;
}

void MakeSpinPlots::DrawIndFit(TString tag, TString mcName){
  TCanvas *cv = new TCanvas(Form("%s_%s",mcName.Data(),tag.Data()));
  
  if(ws->var( Form("Data_%s_INDFIT_%s_Nsig",mcName.Data(),tag.Data()) ) == 0) return;

  RooRealVar* mass = ws->var("mass");
  mass->setBins( (mass->getMax() - mass->getMin())/1.5 ); //enfore 1.5GeV bin width
  RooPlot* frame  = mass->frame();

  tPair lbl(mcName,tag);

  double Ns = ws->var( Form("Data_%s_INDFIT_%s_Nsig",mcName.Data(),tag.Data()) )->getVal();
  double Nb = ws->var( Form("Data_%s_INDFIT_%s_Nbkg",mcName.Data(),tag.Data()) )->getVal();

  double Nblind = ws->data("Data_Combined")->reduce("(mass>100 && mass<119) || (mass>135.5 && mass<170)")->sumEntries(TString("evtcat==evtcat::")+tag);
  double Ntot   = ws->data("Data_Combined")->sumEntries(TString("evtcat==evtcat::")+tag);

  RooFitResult* fitres = (RooFitResult*)ws->obj(Form("Data_%s_INDFIT_fitResult",mcName.Data())); 
  std::cout << fitres << std::endl;
    ws->data("Data_Combined")->reduce(TString("evtcat==evtcat::")+tag)->plotOn(frame,RooFit::LineColor(kWhite),RooFit::MarkerColor(kWhite));
  //Data_Hgg125_INDFIT_EB_0
  ws->pdf(Form("Data_%s_INDFIT_%s",mcName.Data(),tag.Data()))->plotOn(frame, RooFit::FillColor(kGreen),RooFit::VisualizeError(*fitres,2.0));
  ws->pdf(Form("Data_%s_INDFIT_%s",mcName.Data(),tag.Data()))->plotOn(frame, RooFit::FillColor(kYellow),RooFit::VisualizeError(*fitres,1.0));
  ws->pdf(Form("Data_%s_INDFIT_%s",mcName.Data(),tag.Data()))->plotOn(frame, RooFit::LineColor(kRed));
  std::cout << "1" << std::endl;
  ws->pdf(Form("Data_BKGFIT_%s_bkgModel",tag.Data()))->plotOn(frame, RooFit::Normalization(Nb/(Nb+Ns)),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  std::cout << "2" << std::endl;

  ws->data("Data_Combined")->reduce(TString("evtcat==evtcat::")+tag)->plotOn(frame);
  frame->Draw();

  //TLatex *prelim = new TLatex(250,x->GetXmax()-40.,"CMS Preliminary");
  TLatex *prelim = new TLatex(0.12,0.96,"CMS Preliminary");
  TLatex *lum = new TLatex(0.7,0.96,Form("#sqrt{s}=8 TeV  L = %0.1f fb^{-1}",lumi));
  prelim->SetNDC();
  lum->SetNDC();
  prelim->SetTextSize(0.045);
  prelim->SetTextColor(kBlack);
  lum->SetTextSize(0.045);
  lum->SetTextColor(kBlack);

  TLatex *owner = new TLatex(0.6,0.88,"Alex Mott (Nov. 13, 2012)");
  owner->SetNDC();
  owner->SetTextSize(0.045);
  owner->SetTextColor(kBlack);

  TLatex *mu = new TLatex(0.7,0.8,Form("#mu = %0.1f #pm %0.2f", fitMean[lbl].first,fitMean[lbl].second));
  mu->SetNDC();
  mu->SetTextSize(0.045);

  TLatex *sig = new TLatex(0.7,0.72,Form("#sigma_{eff} = %0.1f #pm %0.2f", fitSigEff[lbl].first,fitSigEff[lbl].second));
  sig->SetNDC();
  sig->SetTextSize(0.045);

  float nSig = ws->var( Form("Data_%s_INDFIT_%s_Nsig",mcName.Data(),tag.Data()) )->getVal();
  float nSigErr = ws->var( Form("Data_%s_INDFIT_%s_Nsig",mcName.Data(),tag.Data()) )->getError();

  TLatex *Nsig = new TLatex(0.7,0.64,Form("N_{sig}= %0.1f #pm %0.1f",nSig,nSigErr));
  Nsig->SetNDC();
  Nsig->SetTextSize(0.045);


  frame->addObject(prelim);
  frame->addObject(lum);
  //frame->addObject(owner);
  frame->addObject(mu);
  frame->addObject(sig);
  frame->addObject(Nsig);
  frame->Draw();
  cv->SaveAs( basePath+Form("/mgg-FloatedFraction-%s-%s-%s.png",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/C/mgg-FloatedFraction-%s-%s-%s.C",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/mgg-FloatedFraction-%s-%s-%s.pdf",outputTag.Data(),mcName.Data(),tag.Data()) );
  delete cv;
}

void MakeSpinPlots::DrawSpinBackground(TString tag, TString mcName,bool signal){
  bool drawSM = (smName!="" && smName!=mcName);

  TCanvas cv;
  double thisN  = ws->data(mcName+"_Combined")->reduce(TString("evtcat==evtcat::")+tag)->sumEntries();
  float norm = thisN; //607*lumi/12.*thisN/(totEB+totEE);
  cout << norm <<endl;
  if(signal) norm = ws->data(Form("Data_%s_%s_sigWeight",tag.Data(),mcName.Data()))->sumEntries();
  RooPlot *frame = ws->var("cosT")->frame(0,1,5);

  RooDataSet* bkgWeight = (RooDataSet*)ws->data(Form("Data_%s_%s_bkgWeight",tag.Data(),mcName.Data()));
  RooDataSet* tmp = (RooDataSet*)ws->data("Data_Combined")->reduce(TString("((mass>115 && mass<120) || (mass>130 && mass<135)) && evtcat==evtcat::")+tag);
  tmp->plotOn(frame,RooFit::Rescale(norm/tmp->sumEntries()));
  cout << "b" <<endl;
  ws->pdf(Form("%s_FIT_%s_cosTpdf",mcName.Data(),tag.Data()))->plotOn(frame,RooFit::LineColor(kGreen),RooFit::Normalization(norm/tmp->sumEntries()));
  if(drawSM) ws->pdf(Form("%s_FIT_%s_cosTpdf",smName.Data(),tag.Data()))->plotOn(frame,RooFit::LineColor(kRed),RooFit::Normalization(norm/tmp->sumEntries()));
  cout << "c   " <<bkgWeight <<endl;

  bkgWeight->plotOn(frame,RooFit::Rescale(norm/bkgWeight->sumEntries()),RooFit::MarkerColor(kBlue) );  
  if(signal){
    cout << "d" <<endl;
      
    ws->data(Form("Data_%s_%s_sigWeight",tag.Data(),mcName.Data()))->plotOn(frame,RooFit::MarkerStyle(4));
  }
  cout << "d" <<endl;
  
  frame->SetMaximum(frame->GetMaximum()*(signal?0.8:0.4)*norm/tmp->sumEntries());
  frame->SetMinimum(-1*frame->GetMaximum());
  TLegend l(0.6,0.2,0.95,0.45);
  l.SetFillColor(0);
  l.SetBorderSize(0);
  l.SetHeader(tag);
  l.AddEntry(frame->getObject(0),"Data m#in [115,120]#cup[130,135]","p");
  l.AddEntry(frame->getObject(1),mcName,"l");
  if(drawSM) l.AddEntry(frame->getObject(2),"SM Higgs","l");
  l.AddEntry(frame->getObject(2+drawSM),"background weighted Data","p");
  if(signal) l.AddEntry(frame->getObject(3+drawSM),"signal weighted Data","p");

  cout << "e" <<endl;

  frame->Draw();
  l.Draw("SAME");
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_%s%s_%s_%s.png",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/C/CosThetaDist_%s%s_%s_%s.C",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_%s%s_%s_%s.pdf",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
}

void MakeSpinPlots::DrawSpinSubBackground(TString tag, TString mcName,bool signal){
  bool drawSM = (smName!="" && smName!=mcName);

  TCanvas cv;
  double thisN  = ws->data(mcName+"_Combined")->reduce(TString("evtcat==evtcat::")+tag)->sumEntries();
  float norm = thisN; //607*lumi/12.*thisN/(totEB+totEE);
  tPair lbl(mcName,tag);


  if(signal) norm = nSignal[lbl].first;   //((RooFormulaVar*)ws->obj(Form("Data_%s_INDFIT_%s_Nsig",mcName.Data(),tag.Data())) )->getVal();
  RooPlot *frame = ws->var("cosT")->frame(0,1,10);

  RooDataSet* tmp = (RooDataSet*)ws->data("Data_Combined")->reduce(TString("((mass>115 && mass<120) || (mass>130 && mass<135)) && evtcat==evtcat::")+tag);
  tmp->plotOn(frame,RooFit::Rescale(norm/tmp->sumEntries()));

  ws->pdf(Form("%s_FIT_%s_cosTpdf",mcName.Data(),tag.Data()))->plotOn(frame,RooFit::LineColor(kGreen),RooFit::Normalization(norm/tmp->sumEntries()));
  if(drawSM) ws->pdf(Form("%s_FIT_%s_cosTpdf",smName.Data(),tag.Data()))->plotOn(frame,RooFit::LineColor(kRed),RooFit::Normalization(norm/tmp->sumEntries()));
  if(signal){
    RooDataHist *h = (RooDataHist*)ws->data( Form("Data_%s_%s_bkgSub_cosT",mcName.Data(),tag.Data()) );
    h->plotOn(frame,RooFit::MarkerStyle(4));
    std::cout << "Nsig: " << h->sumEntries() << std::endl;
  }
  
  frame->SetMaximum(frame->GetMaximum()*(signal?3.:1.2)*norm/tmp->sumEntries());
  frame->SetMinimum(-1*frame->GetMaximum());
  TLegend l(0.6,0.2,0.95,0.45);
  l.SetFillColor(0);
  l.SetBorderSize(0);
  l.SetHeader(tag);
  l.AddEntry(frame->getObject(0),"Data m#in [115,120]#cup[130,135]","p");
  l.AddEntry(frame->getObject(1),mcName,"l");
  if(drawSM) l.AddEntry(frame->getObject(2),"SM Higgs","l");
  if(signal) l.AddEntry(frame->getObject(2+drawSM),"bkg-subtracted Data","p");
  
  frame->Draw();
  l.Draw("SAME");
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s_%s.png",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/C/CosThetaDist_SimpleSub_%s%s_%s_%s.C",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s_%s.pdf",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
}

void MakeSpinPlots::DrawSpinSubTotBackground(TString mcName,bool signal){
  bool drawSM = (smName!="" && smName!=mcName);

  TCanvas cv;
  double thisN  = ws->data(mcName+"_Combined")->sumEntries();
  float norm = thisN;


  if(signal) norm = ws->var(Form("Data_%s_FULLFIT_Nsig",mcName.Data()))->getVal();
  RooPlot *frame = ws->var("cosT")->frame(0,1,10);

  RooDataSet* tmp = (RooDataSet*)ws->data(Form("Data_Combined"))->reduce("(mass>115 && mass<120) || (mass>130 && mass<135)");
  tmp->plotOn(frame,RooFit::Rescale(norm/tmp->sumEntries()));

  ws->pdf(Form("%s_FIT_cosTpdf",mcName.Data()))->plotOn(frame,RooFit::LineColor(kGreen),RooFit::Normalization(norm/tmp->sumEntries()));
  if(drawSM)  ws->pdf(Form("%s_FIT_cosTpdf",smName.Data()))->plotOn(frame,RooFit::LineColor(kRed),RooFit::Normalization(norm/tmp->sumEntries()));
  if(signal){
    RooDataHist *h = (RooDataHist*)ws->data( Form("Data_%s_Combined_bkgSub_cosT",mcName.Data()) );
    h->plotOn(frame,RooFit::MarkerStyle(4));
    std::cout << "Nsig: " << h->sumEntries() << std::endl;
  }

  
  frame->SetMaximum(frame->GetMaximum()*(signal?2.:1.2)*norm/tmp->sumEntries());
  frame->SetMinimum(-1*frame->GetMaximum());
  TLegend l(0.6,0.2,0.95,0.45);
  l.SetFillColor(0);
  l.SetBorderSize(0);
  l.SetHeader("Combined");
  l.AddEntry(frame->getObject(0),"Data m#in [115,120]#cup[130,135]","p");
  l.AddEntry(frame->getObject(1),mcName,"l");
  if(drawSM) l.AddEntry(frame->getObject(2),"SM Higgs","l");
  if(signal) l.AddEntry(frame->getObject(2+drawSM),"bkg-subtracted Data","p");
  
  frame->Draw();
  l.Draw("SAME");
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s.png",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/C/CosThetaDist_SimpleSub_%s%s_%s.C",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s.pdf",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data()) );
}


void MakeSpinPlots::PlotSignalFits(TString tag, TString mcName,TString cosThetaBin){
  TCanvas cv;
  TString cat=tag;
  if(cosThetaBin!="") tag = tag+"_"+cosThetaBin;

  float mean = ws->var(Form("%s_FIT_%s_mean",mcName.Data(),tag.Data()))->getVal();
  RooPlot *frame = ws->var("mass")->frame(105,140,70);//mean-10,mean+10,40);
  RooAbsData *d = ws->data(mcName+"_Combined")->reduce(TString("evtcat==evtcat::")+cat);
  if(cosThetaBin!=""){
    TObjArray *arr = cosThetaBin.Tokenize("_");
    float low  = atof(arr->At(1)->GetName());
    float high = atof(arr->At(2)->GetName());
    d = d->reduce( Form("cosT < %0.2f && cosT >= %0.2f",high,low) );
    delete arr;
  }

  d->plotOn(frame);
  RooFitResult *res = (RooFitResult*)ws->obj(Form("%s_FIT_%s_fitResult",mcName.Data(),tag.Data()));
  RooAbsPdf * pdf = ws->pdf(Form("%s_FIT_%s",mcName.Data(),tag.Data())); //signal model
  std::cout << pdf << "\t" << res << std::endl;
  pdf->plotOn(frame,RooFit::FillColor(kGreen),RooFit::VisualizeError(*res,2.0));
  pdf->plotOn(frame,RooFit::FillColor(kYellow),RooFit::VisualizeError(*res,1.0));
  pdf->plotOn(frame,RooFit::LineColor(kRed));
  d->plotOn(frame); //data
  
  tPair lbl(mcName,tag);

  TLatex *prelim = new TLatex(0.18,0.9,"CMS Preliminary Simulation");
  TLatex *sigL  = new TLatex(0.18,0.6,Form("#sigma_{eff} = %0.2f GeV",fitSigEff[lbl].first,fitSigEff[lbl].second));
  prelim->SetNDC();
  sigL->SetNDC();
  prelim->SetTextSize(0.05);
  sigL->SetTextSize(0.05);
  
  frame->addObject(prelim);
  frame->addObject(sigL);
  frame->Draw();
  cv.SaveAs(basePath+Form("/signalModels/sig_%s_%s_%s.png",mcName.Data(),outputTag.Data(),tag.Data()));
  cv.SaveAs(basePath+Form("/signalModels/C/sig_%s_%s_%s.C",mcName.Data(),outputTag.Data(),tag.Data()));
  cv.SaveAs(basePath+Form("/signalModels/sig_%s_%s_%s.pdf",mcName.Data(),outputTag.Data(),tag.Data()));

}

void MakeSpinPlots::setStyle(){

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
  vecbosStyle->SetPadRightMargin(0.05);
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

void MakeSpinPlots::printAll(){
  std::vector<TString>::const_iterator mcIt = mcNames.begin();
  for(; mcIt != mcNames.end(); mcIt++){
    std::cout << "\n\n" << *mcIt << std::endl;
    printYields(*mcIt);
  }

}
void MakeSpinPlots::printYields(const char* mcType){
  RooRealVar * tot = ws->var(Form("Data_%s_FULLFIT_Nsig",mcType));
  if(tot==0) return;

  cout << "Total Yield:  " << tot->getVal() << "  +-  " << tot->getError() <<endl;
  cout << "Category Yields: CONSTRAINED FIT " << endl;
  for(int i=0;i<catNames.size();i++){
    RooRealVar *f = ws->var( Form("Data_%s_FULLFIT_%s_fsig",mcType, catNames.at(i).Data()) );
    cout << "\t" << catNames.at(i) <<":   " << tot->getVal()*f->getVal() << "  +-  " << tot->getError()*f->getVal() <<endl;
  }
  cout << "\nCategory Yields: INDEPENDENT FIT " << endl;
  for(int i=0;i<catNames.size();i++){
    RooRealVar *ind = ws->var( Form("Data_%s_INDFIT_%s_Nsig",mcType, catNames.at(i).Data()) );
    if(ind==0) continue;
    cout << "\t" << catNames.at(i) <<":   " << ind->getVal() << "  +-  " << ind->getError() <<endl;
  }

  float exp = ws->data(Form("%s_Combined",mcType))->sumEntries();///total * 607*lumi/12.;
  cout << endl << "Expected Events:  "  << exp << endl;
  cout << "Expected Yields Per Category: " <<endl;
  for(int i=0;i<catNames.size();i++){ 
    RooRealVar *f = ws->var( Form("Data_%s_FULLFIT_%s_fsig",mcType, catNames.at(i).Data()) );
    cout << "\t" << catNames.at(i) <<":   " << exp*f->getVal() <<endl;
  }

  cout << "mu:  " << tot->getVal()/exp << "  +-  "
       << tot->getError()/exp <<endl;
  
  for(int i=0;i<catNames.size();i++){
    RooRealVar *ind = ws->var( Form("Data_%s_INDFIT_%s_Nsig",mcType, catNames.at(i).Data()) );
    if(ind==0) continue;
    RooRealVar *f = ws->var( Form("Data_%s_FULLFIT_%s_fsig",mcType, catNames.at(i).Data()) );
    cout << "\t" << catNames.at(i) <<":   " << ind->getVal()/(exp*f->getVal()) << "  +-  " << ind->getError()/(exp*f->getVal()) <<endl;
  }
  MakeChannelComp(mcType);
}

void MakeSpinPlots::MakeChannelComp(const char* mcType){
  TGraphErrors graph(catNames.size());

  RooRealVar mu("mu","",-50,50);

  float exp = ws->data(Form("%s_Combined",mcType))->sumEntries();///total * 607*lumi/12.;

  TH1F frame("frame","",catNames.size(),0,catNames.size());

  //graph.GetXaxis()->SetNdivisions(catNames.size());

  float min=99999,max=-99999;
  for(int i=0;i<catNames.size();i++){
    RooRealVar *ind = ws->var( Form("Data_%s_INDFIT_%s_Nsig",mcType, catNames.at(i).Data()) );
    RooRealVar *f = ws->var( Form("Data_%s_FULLFIT_%s_fsig",mcType, catNames.at(i).Data()) );
    float mu = ind->getVal()/exp/f->getVal();
    float muE = ind->getError()/exp/f->getVal();
    graph.SetPoint(i,i+0.5,mu);
    graph.SetPointError(i,0,muE);

    if(mu-muE < min) min = mu-muE;
    if(mu+muE > max) max = mu+muE;
    
    //graph.GetXaxis()->SetBinLabel(i+1,catNames.at(i));
  }

  TF1 fit("fit","[0]",0+0.5,catNames.size()+0.5);
  graph.Fit(&fit,"MNE");
  float mean  = fit.GetParameter(0);
  float meanE = fit.GetParError(0);
  
  frame.SetAxisRange(min-1.5,max+0.5,"Y");
  frame.SetYTitle("Fitted #sigma/#sigma_{SM}");
  frame.SetXTitle("Category");
  TCanvas cv;
  frame.Draw();

  TBox err(0,mean-meanE,catNames.size(),mean+meanE);
  err.SetFillColor(kGreen);

  frame.Draw();
  err.Draw("SAME");
  TLine mLine(0,mean,catNames.size(),mean);
  mLine.Draw("SAME");
  graph.Draw("PSAME");
  

  cv.SaveAs(basePath+Form("/ChannelComp_%s_%s.png",mcType,outputTag.Data()));
  cv.SaveAs(basePath+Form("/ChannelComp_%s_%s.pdf",mcType,outputTag.Data()));
}

void MakeSpinPlots::setupDir(){
  mkdir(basePath.Data(),S_IRWXU|S_IRGRP|S_IXGRP);

  basePath = basePath+"/"+outputTag+"/";
  mkdir(basePath.Data(),S_IRWXU|S_IRGRP|S_IXGRP);    
  TString path = basePath;
  //C-files
  path = basePath+"C";
  mkdir((path.Data()),S_IRWXU|S_IRGRP|S_IXGRP);

  //signal models
  path = basePath+"signalModels";
  mkdir((path.Data()),S_IRWXU|S_IRGRP|S_IXGRP);
  path = basePath+"signalModels/C/";
  mkdir((path.Data()),S_IRWXU|S_IRGRP|S_IXGRP);


  // cos(theta)
  path = basePath+"cosThetaPlots";
  mkdir((path.Data()),S_IRWXU|S_IRGRP|S_IXGRP);
  path = basePath+"cosThetaPlots/C/";
  mkdir((path.Data()),S_IRWXU|S_IRGRP|S_IXGRP);

  isSetup = true;
}
