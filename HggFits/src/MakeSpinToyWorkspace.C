#include "MakeSpinToyWorkspace.h"

MakeSpinToyWorkspace::MakeSpinToyWorkspace(TString fileName,TString wsName):MakeSpinSeparationTest(fileName,wsName){
  mass = ws->var("mass");
  saveWorkspaces=false;
}

void MakeSpinToyWorkspace::save(TString filename){
  TFile *f = new TFile(filename,"RECREATE");
  this->makeForCombineTool()->Write();
  f->Close();
}

void MakeSpinToyWorkspace::setup(TString t,float l){
  tag=t;
  lumi=l;
  hggMassPdf = ws->pdf( Form("Hgg125_FIT_%s",tag.Data()) );
  rsgMassPdf = ws->pdf( Form("RSG125_FIT_%s",tag.Data()) );
  bkgMassPdf = ws->pdf( Form("Data_Hgg125_BKGFIT_%s_bkgModel",tag.Data()) );


  hggPdf = (RooHistPdf*)ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",tag.Data()));
  rsgPdf = (RooHistPdf*)ws->pdf(Form("RSG125_FIT_%s_cosTpdf",tag.Data()));
  getBackgroundPdf(tag);

  Nbkg = ws->var(Form("Data_Hgg125_BKGFIT_%s_Nbkg",tag.Data()))->getVal();
  Nsig = (int)getExpEvents(lumi,tag);

}

void MakeSpinToyWorkspace::generateN(int Ntoy){
  for(int i=0;i<Ntoy;i++){
    generateToyWorkspace(true);
    generateToyWorkspace(false);
  }
}



void MakeSpinToyWorkspace::generateToyWorkspace(bool doRSG){
  //RooProdPdf bkgGen("bkgGen","@0*@1",RooArgList(*bkgMassPdf,*bkgGenPdf));
  //RooDataSet *M = bkgGen.generate(RooArgSet(*mass,*cosT),Nbkg,RooFit::Extended());
  TRandom3 rng(0);
  int thisNbkg = (int)rng.Poisson(Nbkg);
  int thisNsig = (int)rng.Poisson(Nsig);

  RooDataSet *M = bkgMassPdf->generate(*mass,thisNbkg);
  RooDataSet *C = bkgGenPdf->generate(*cosT,thisNbkg);

  //RooProdPdf * sigGen;
  RooDataSet *mcM, *mcC;
  if(doRSG){
    //sigGen = new RooProdPdf("sigGen","@0*@1",RooArgList(*rsgMassPdf,*rsgGenPdf));
    mcM = rsgMassPdf->generate(*mass,thisNsig);
    mcC = rsgGenPdf->generate(*cosT,thisNsig);
  }else{
    //sigGen = new RooProdPdf("sigGen","@0*@1",RooArgList(*hggMassPdf,*hggGenPdf));
    mcM = hggMassPdf->generate(*mass,thisNsig);
    mcC = hggGenPdf->generate(*cosT,thisNsig);
  }
  //RooDataSet *mcM = sigGen->generate(RooArgSet(*mass,*cosT),Nsig,RooFit::Extended());

  mcM->SetName(Form("ToySigData_%s",tag.Data()));
  M->merge(C);
  mcM->merge(mcC);
  M->append(*mcM);

  M->SetName(Form("Data_%s",tag.Data()));

  //delete sigGen;

  RooWorkspace *toyws = new RooWorkspace("ws","");
  toyws->import(*M);

  toyws->import( *hggPdf );
  toyws->import( *rsgPdf );
  toyws->import( *hggMassPdf );
  toyws->import( *rsgMassPdf );

  //ws->import( *mcM );

  MakeSpinFits fits("","");

  fits.setWorkspace(toyws);  
  fits.MakeBackgroundFit("Hgg125",tag,125,2,false);
  fits.MakeBackgroundFit("RSG125",tag,125,2,false);

  fits.MakeBackgroundFitCosTBin("Hgg125",tag,-1,-0.3);
  fits.MakeBackgroundFitCosTBin("RSG125",tag,-1,-0.3);

  fits.MakeBackgroundFitCosTBin("Hgg125",tag,-0.3,0.3);
  fits.MakeBackgroundFitCosTBin("RSG125",tag,-0.3,0.3);

  fits.MakeBackgroundFitCosTBin("Hgg125",tag,0.3,1);
  fits.MakeBackgroundFitCosTBin("RSG125",tag,0.3,1);

  float N[3],Ne[3];
  N[0] = fits.getWorkspace()->var( Form("Data_Hgg125_FIT_-1.0_cosT_-0.3_%s_Nsig",tag.Data()) )->getVal();
  N[1] = fits.getWorkspace()->var( Form("Data_Hgg125_FIT_-0.3_cosT_0.3_%s_Nsig",tag.Data()) )->getVal();
  N[2] = fits.getWorkspace()->var( Form("Data_Hgg125_FIT_0.3_cosT_1.0_%s_Nsig",tag.Data()) )->getVal();

  Ne[0] = fits.getWorkspace()->var( Form("Data_Hgg125_FIT_-1.0_cosT_-0.3_%s_Nsig",tag.Data()) )->getError();
  Ne[1] = fits.getWorkspace()->var( Form("Data_Hgg125_FIT_-0.3_cosT_0.3_%s_Nsig",tag.Data()) )->getError();
  Ne[2] = fits.getWorkspace()->var( Form("Data_Hgg125_FIT_0.3_cosT_1.0_%s_Nsig",tag.Data()) )->getError();

  float expHgg[3];
  expHgg[0] = getExpEventsCosT(lumi,"Hgg125",tag,-1.0,-0.3);
  expHgg[1] = getExpEventsCosT(lumi,"Hgg125",tag,-0.3, 0.3);
  expHgg[2] = getExpEventsCosT(lumi,"Hgg125",tag, 0.3, 1.0);
  
  float expRSG[3];
  expRSG[0] = getExpEventsCosT(lumi,"RSG125",tag,-1.0,-0.3);
  expRSG[1] = getExpEventsCosT(lumi,"RSG125",tag,-0.3, 0.3);
  expRSG[2] = getExpEventsCosT(lumi,"RSG125",tag, 0.3, 1.0);
  
  double lHgg = TMath::Log(prob(expHgg,N,Ne,3));
  double lRSG = TMath::Log(prob(expRSG,N,Ne,3));

  if(doRSG==false){
    nll->setVal(lHgg);
    hgg_ds_Thgg->add(RooArgSet(*nll));
    nll->setVal(lRSG);
    rsg_ds_Thgg->add(RooArgSet(*nll));
    S->setVal(2*(lRSG-lHgg) );
    s_ds_Thgg->add(RooArgSet(*S));
  }else{
    nll->setVal(lHgg);
    hgg_ds_Trsg->add(RooArgSet(*nll));
    nll->setVal(lRSG);
    rsg_ds_Trsg->add(RooArgSet(*nll));
    S->setVal(2*(lRSG-lHgg) );
    s_ds_Trsg->add(RooArgSet(*S));
  }                             
    
  if(saveWorkspaces) toyWorkspaces.Add(toyws);
  else toyws->Delete();
}

float MakeSpinToyWorkspace::getExpEventsCosT(float lumi, TString MC, TString cat,float minCosT, float maxCosT){
  double totEB  = ws->var(Form("%s_EB_totalEvents",MC.Data()))->getVal();
  double totEE  = ws->var(Form("%s_EE_totalEvents",MC.Data()))->getVal();

  double thisN = ws->data(Form("%s_%s",MC.Data(),cat.Data()))->reduce(Form("cosT>=%f && cosT<%f",minCosT,maxCosT))->sumEntries();

  return thisN/(totEB+totEE)*lumi/12*607; //607 events in 12/fb @ 8 TeV
}

double MakeSpinToyWorkspace::prob(float *exp, float *obs, float *err,int N){
  double chi2=0;
  int ndof=0;
  for(int i=0;i<N-1;i++){
    double expRatio = exp[i+1]/exp[i];
    double obsRatio = obs[i+1]/obs[i];
    double obsRatioE = obsRatio*TMath::Sqrt( err[i+1]*err[i+1]/obs[i+1]/obs[i+1]+err[i]*err[i]/obs[i]/obs[i]);
    if(obsRatioE){
      chi2 += TMath::Power( (expRatio-obsRatio)/obsRatioE, 2);
      ndof++;
    }
  }

  return TMath::Prob(chi2,ndof);
}
