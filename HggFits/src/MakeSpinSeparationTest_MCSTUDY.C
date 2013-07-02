#include "MakeSpinWorkspace.C"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"

#include <vector>

class MakeSpinSeparationTest{
public:
  MakeSpinSeparationTest(TString fileName,TString wsName="cms_hgg_spin_workspace");
  RooMCStudy* runMCStudy(int Ntoys,float lumi,TString cat,bool doHgg);

  float getExpEvents(float lumi, TString cat);

  std::vector<TH1F*> hggHists,rsgHists;

  RooWorkspace *ws;
  RooRealVar* cosT;
};

MakeSpinSeparationTest::MakeSpinSeparationTest(TString fileName,TString wsName){
  TFile *f = new TFile(fileName);
  ws = (RooWorkspace*)f->Get(wsName);
  cosT = ws->var("cosT");
  cosT->setBins(5);
}

RooMCStudy* MakeSpinSeparationTest::runMCStudy(int Ntoys,float lumi,TString cat,bool doHgg){
  RooDataSet* hggDS = (RooDataSet*)ws->data(Form("Hgg125_%s",cat.Data()));
  RooDataSet* rsgDS = (RooDataSet*)ws->data(Form("RSG125_%s",cat.Data()));

  RooKeysPdf *hggPdf = new RooKeysPdf(Form("Hgg125_%s_KDE",cat.Data()),"",*cosT,*hggDS);
  RooKeysPdf *rsgPdf = new RooKeysPdf(Form("RSG125_%s_KDE",cat.Data()),"",*cosT,*rsgDS);

  cout << cosT << "  " << hggPdf << "  " << rsgPdf <<endl;

  RooChi2MCSModule chi2mod;


  RooMCStudy *study;

  if(doHgg) study = new RooMCStudy(*hggPdf,*cosT,RooFit::Binned());
  else      study = new RooMCStudy(*hggPdf,*cosT,RooFit::FitModel(*rsgPdf),RooFit::Binned());

  //study->addModule(chi2mod);

  int N = (int)getExpEvents(lumi,cat);

  study->generateAndFit(Ntoys,N,kTRUE);
  return study;

}

float MakeSpinSeparationTest::getExpEvents(float lumi, TString cat){
  double totEB  = ws->var("Hgg125_EB_totalEvents")->getVal();
  double totEE  = ws->var("Hgg125_EE_totalEvents")->getVal();

  double thisN  = ws->data(Form("Hgg125_%s",cat.Data()))->sumEntries();
  
  return thisN/(totEB+totEE)*lumi/12*607; //607 events in 12/fb @ 8 TeV
}
