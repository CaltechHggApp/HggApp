#include "MakeSpinSeparationTest.h"
#include "MakeSpinWorkspace.C"


TTree* MakeSpinSeparationTest::makeForCombineTool(){
  TTree * out = new TTree("q","");
  float q;
  int type;

  out->Branch("q",&q);
  out->Branch("type",&type,"type/I");

  Long64_t iEntry=-1;
  type=1; //SM HIggs
  const RooArgSet *set;
  while( (set=s_ds_Thgg->get(++iEntry)) ){ 
    q = ((RooRealVar*)set->find("S"))->getVal();
    out->Fill();
  }
  iEntry=-1;
  type=-1;
  while( (set=s_ds_Trsg->get(++iEntry)) ){ 
    q = ((RooRealVar*)set->find("S"))->getVal();
    out->Fill();
  }
  return out;
}

void MakeSpinSeparationTest::getBackgroundPdf(TString cat){
  RooDataSet *tmp  = (RooDataSet*)(ws->data(Form("Data_%s",cat.Data()))->reduce("(mass>110 && mass<120) || (mass>130 && mass<140)"));
  RooDataSet *tmpH = (RooDataSet*)ws->data(Form("Hgg125_%s",cat.Data()));
  RooDataSet *tmpR = (RooDataSet*)ws->data(Form("RSG125_%s",cat.Data()));
  bkgGenPdf = new RooKeysPdf("bkgGenPdf","",*cosT,*tmp);
  hggGenPdf = new RooKeysPdf("hggGenPdf","",*cosT,*tmpH);
  rsgGenPdf = new RooKeysPdf("rsgGenPdf","",*cosT,*tmpR);

  RooDataHist *tmpHist = new RooDataHist("tmpHist","",RooArgSet(*cosT),*tmp);
  bkgPdf = new RooHistPdf("bkgPdf","",RooArgSet(*cosT),*tmpHist);
}

MakeSpinSeparationTest::MakeSpinSeparationTest(TString fileName,TString wsName){
  TFile *f = new TFile(fileName);
  ws = (RooWorkspace*)f->Get(wsName);
  cosT = ws->var("cosT");
  cosT->setBins(5);
  nll = new RooRealVar("nll","",0,0,1e5);
  S = new RooRealVar("S","",0,-100,100);

  hgg_ds_Thgg = new RooDataSet("Hgg125_nll_Thgg","",RooArgSet(*nll));
  rsg_ds_Thgg = new RooDataSet("RSG125_nll_Thgg","",RooArgSet(*nll));
  s_ds_Thgg   = new RooDataSet("S_nll_Thgg","",RooArgSet(*S));

  hgg_ds_Trsg = new RooDataSet("Hgg125_nll_Trsg","",RooArgSet(*nll));
  rsg_ds_Trsg = new RooDataSet("RSG125_nll_Trsg","",RooArgSet(*nll));
  s_ds_Trsg   = new RooDataSet("S_nll_Trsg","",RooArgSet(*S));
}



void MakeSpinSeparationTest::runMCStudy(int Ntoys,float lumi,TString cat){
  RooDataSet* hggDS = (RooDataSet*)ws->data(Form("Hgg125_%s",cat.Data()));
  RooDataSet* rsgDS = (RooDataSet*)ws->data(Form("RSG125_%s",cat.Data()));

  hggPdf = (RooHistPdf*)ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",cat.Data()));
  rsgPdf = (RooHistPdf*)ws->pdf(Form("RSG125_FIT_%s_cosTpdf",cat.Data()));
  getBackgroundPdf(cat);
  //  RooKeysPdf *hggPdf = new RooKeysPdf(Form("Hgg125_%s_KDE",cat.Data()),"",*cosT,*hggDS);
  //  RooKeysPdf *rsgPdf = new RooKeysPdf(Form("RSG125_%s_KDE",cat.Data()),"",*cosT,*rsgDS);

  cout << cosT << "  " << hggPdf << "  " << rsgPdf <<endl;


  N = (int)getExpEvents(lumi,cat);
  float bkgE = sqrt(getBkgEvents(lumi,cat));
  for(int iToy=0;iToy<Ntoys;iToy++){
    std::pair<double,double> n1 = getNLL(N,true,bkgE);
    std::pair<double,double> n2 = getNLL(N,false,bkgE);
    nll->setVal(n1.first);
    hgg_ds_Thgg->add(RooArgSet(*nll));
    nll->setVal(n1.second);
    rsg_ds_Thgg->add(RooArgSet(*nll));
    S->setVal(2*(n1.second-n1.first));
    s_ds_Thgg->add(RooArgSet(*S));

    nll->setVal(n2.first);
    hgg_ds_Trsg->add(RooArgSet(*nll));
    nll->setVal(n2.second);
    rsg_ds_Trsg->add(RooArgSet(*nll));
    S->setVal(2*(n2.first-n2.second));
    s_ds_Trsg->add(RooArgSet(*S));
  }

  RooDataSet* data = (RooDataSet*)ws->data("Hgg125_sigWeight_EB_0_0");
  std::pair<double,double> Sdata = getDataNLL(hggPdf,rsgPdf,cosT,data);
  std::cout << "Data Value:  " << 2*(Sdata.second-Sdata.first) <<endl;

}

float MakeSpinSeparationTest::getExpEvents(float lumi, TString cat){
  double totEB  = ws->var("Hgg125_EB_totalEvents")->getVal();
  double totEE  = ws->var("Hgg125_EE_totalEvents")->getVal();

  double thisN  = ws->data(Form("Hgg125_%s",cat.Data()))->sumEntries();
  
  return thisN/(totEB+totEE)*lumi/12*607; //607 events in 12/fb @ 8 TeV
}

float MakeSpinSeparationTest::getBkgEvents(float lumi, TString cat){
  RooAbsPdf * bkg = ws->pdf(Form("Data_Hgg125_BKGFIT_%s_bkgModel",cat.Data()));

  double s1   = ws->var(Form("Data_Hgg125_FIT_%s_sigma1",cat.Data()))->getVal();
  double s2   = ws->var(Form("Data_Hgg125_FIT_%s_sigma2",cat.Data()))->getVal();
  double f    = ws->var(Form("Data_Hgg125_FIT_%s_f",cat.Data()))->getVal();
  double se = f*s1+(1-f)*s2;
  double Nbkg  = ws->var(Form("Data_Hgg125_FIT_%s_Nbkg",cat.Data()))->getVal();

  RooRealVar range("range","",125-se,125+se);
  RooRealVar all("all","",100,180);
  float BkgInRange = bkg->createIntegral(range)->getVal()/bkg->createIntegral(all)->getVal()*Nbkg;
 
  return BkgInRange*lumi/11.2;
}

std::pair<double,double> MakeSpinSeparationTest::getNLL(float Nev,bool doHgg,float NbkgErr){
  TRandom3 rng(0);
  RooDataSet* ds;
  if(doHgg) ds = hggGenPdf->generate(*cosT,Nev,RooFit::Extended());
  else ds = rsgGenPdf->generate(*cosT,Nev,RooFit::Extended());
  if(!bkgPdf){ //background free version
    hggPdf->fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0));
    RooFitResult* res1 = hggPdf->fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2));
    rsgPdf->fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0));
    RooFitResult* res2 = rsgPdf->fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2));
    return std::pair<double,double>(res1->minNll(),res2->minNll());
  }  else {
    std::cout << "Generating bkg dataset" <<std::endl;
    int Ngen = fabs((int)rng.Gaus(0,NbkgErr));
    if(Ngen){
      std::cout << Ngen <<std::endl;
      RooDataSet* bkg = bkgGenPdf->generate(*cosT,Ngen);
      ds->append(*bkg);
    }

    RooRealVar NbkgH("NbkgH","",0,0,4*NbkgErr);
    RooRealVar NsigH("NsigH","",0,0,1e8);
    RooRealVar NbkgR("NbkgR","",0,0,4*NbkgErr);
    RooRealVar NsigR("NsigR","",0,0,1e8);

    RooAddPdf pH("pH","",RooArgList(*hggPdf,*bkgPdf),RooArgList(NsigH,NbkgH));
    RooAddPdf pR("pR","",RooArgList(*rsgPdf,*bkgPdf),RooArgList(NsigR,NbkgR));

    RooRealVar bErr("bErr","",(NbkgErr<1?1:NbkgErr));
    RooRealVar bMean("bMean","",0.);
    RooGaussian bkgConstH("bkgConstH","",NbkgH,bMean,bErr);
    RooGaussian bkgConstR("bkgConstR","",NbkgR,bMean,bErr);

    pH.fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0),RooFit::ExternalConstraints(bkgConstH));
    RooFitResult* res1 = pH.fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2),RooFit::ExternalConstraints(bkgConstH));
    res1->Print("V");
    pR.fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0),RooFit::ExternalConstraints(bkgConstR));
    RooFitResult* res2 = pR.fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2),RooFit::ExternalConstraints(bkgConstR));
    return std::pair<double,double>(res1->minNll(),res2->minNll());    
  }

}

std::pair<double,double> MakeSpinSeparationTest::getDataNLL(RooAbsPdf* hggPdf, RooAbsPdf* rsgPdf,RooRealVar* var,RooAbsData* ds,float NbkgErr,float Nsig,float NsigErr){
  if(NbkgErr==0){
    hggPdf->fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0));
    RooFitResult* res1 = hggPdf->fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2));
    rsgPdf->fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0));
    RooFitResult* res2 = rsgPdf->fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2));
    return std::pair<double,double>(res1->minNll(),res2->minNll());
  }


  RooRealVar NbkgH("NbkgH","",0,0,0*NbkgErr);
  RooRealVar NsigH("NsigH","",Nsig,0,1e8);
  RooRealVar NbkgR("NbkgR","",0,0,0*NbkgErr);
  RooRealVar NsigR("NsigR","",Nsig,0,1e8);
  
  RooAddPdf pH("pH","",RooArgList(*hggPdf,*bkgPdf),RooArgList(NsigH,NbkgH));
  RooAddPdf pR("pR","",RooArgList(*rsgPdf,*bkgPdf),RooArgList(NsigR,NbkgR));
  
  RooRealVar bErr("bErr","",(NbkgErr<1?1:NbkgErr));
  RooRealVar bMean("bMean","",0.);
  RooGaussian bkgConstH("bkgConstH","",NbkgH,bMean,bErr);
  RooGaussian bkgConstR("bkgConstR","",NbkgR,bMean,bErr);

  RooRealVar sErr("sErr","",NsigErr);
  RooRealVar sMean("sMean","",Nsig);
  
  RooGaussian sigConstH("sigConstH","",NsigH,sMean,sErr);
  RooGaussian sigConstR("sigConstR","",NsigR,sMean,sErr);

  pH.fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0),RooFit::ExternalConstraints(RooArgList(bkgConstH,sigConstH)));
  RooFitResult* res1 = pH.fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2),RooFit::ExternalConstraints(RooArgList(bkgConstH,sigConstH)));
  pR.fitTo(*ds,RooFit::Extended(),RooFit::Strategy(0),RooFit::ExternalConstraints(RooArgList(bkgConstR,sigConstR)));
  RooFitResult* res2 = pR.fitTo(*ds,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2),RooFit::ExternalConstraints(RooArgList(bkgConstR,sigConstR)));
  return std::pair<double,double>(res1->minNll(),res2->minNll());
}
