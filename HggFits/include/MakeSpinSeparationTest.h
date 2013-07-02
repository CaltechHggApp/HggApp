#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooGlobalFunc.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include "RooStats/SPlot.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooBernstein.h"
#include "RooLinkedListIter.h"
#include "RooCBShape.h"


#include <vector>
#include <TRandom3.h>

class MakeSpinSeparationTest{
public:
  MakeSpinSeparationTest(TString fileName,TString wsName="cms_hgg_spin_workspace");
  void runMCStudy(int Ntoys,float lumi,TString cat);


  float getExpEvents(float lumi, TString cat);
  float getBkgEvents(float lumi, TString cat);

  std::pair<double,double> getNLL(float Nev,bool genHgg,float NbkgErr=0);

  RooWorkspace *ws;
  RooRealVar* cosT;
  RooRealVar *nll,*S;
  RooDataSet *hgg_ds_Thgg,*rsg_ds_Thgg,*s_ds_Thgg;
  RooDataSet *hgg_ds_Trsg,*rsg_ds_Trsg,*s_ds_Trsg;

  RooHistPdf *rsgPdf,*hggPdf,*bkgPdf;
  RooKeysPdf *bkgGenPdf,*hggGenPdf,*rsgGenPdf;
  int N;

  std::pair<double,double> getDataNLL(RooAbsPdf* hggPdf, RooAbsPdf* rsgPdf,RooRealVar* var,RooAbsData* ds,float NbkgErr=0,float Nsig=0,float NsigErr=0);

  void getBackgroundPdf(TString cat);
  TTree* makeForCombineTool();
};
