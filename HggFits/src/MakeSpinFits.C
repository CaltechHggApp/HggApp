#include "MakeSpinFits.h"
#include "subtract.cc"

#include "RooTrace.h"
#define NUM_CPU 1

#define __DO_TRACE 0
using namespace std;
//default constructor
MakeSpinFits::MakeSpinFits():
  addSWeight(true),
  useCB(false),
  use2DIntegralMorph(false),
  cosTbinEdges(0)
{
  ws=0;
  setBkgFit(MakeSpinFits::kExp);
  setMeanRange(124,126,125);
}

MakeSpinFits::MakeSpinFits(TString inputFileName, TString outputFileName):
  addSWeight(true),
  useCB(false),
  use2DIntegralMorph(false)
{
  if(inputFileName != ""){
    //opens the input file and gets the input workspace
    inputFile = new TFile(inputFileName);
    ws = ((RooWorkspace*)inputFile->Get("cms_hgg_spin_workspace"));
    ws->Print();
    //extract MC labels from the input workspace
    getLabels("labels",&mcLabel,ws);
    //extract category labels from the input workspace
    getLabels("evtcat",&catLabels,ws);
  }
  if(outputFileName != ""){
    //opens the output file
    outputFile = new TFile(outputFileName,"RECREATE");
    outputFile->cd();
    ws->Write(ws->GetName(),TObject::kWriteDelete);
  }
  //default fit type
  setBkgFit(MakeSpinFits::kExp);
  setMeanRange(123,126,125);

  //default cos(theta) binning
  float tmp[5] = {0, 0.2, 0.4, 0.6, 1};
  setCosTBins(5,tmp);
}

MakeSpinFits::~MakeSpinFits(){
  RooTrace::dump(cout);
}
//get the labels (MC names or category names) from the workspace and put them in lblVec
void MakeSpinFits::getLabels(const char *varName, std::vector<TString> *lblVec,RooWorkspace* w){
  RooCategory* labels = ((RooCategory*)w->obj(varName));
  lblVec->clear();
  if(labels==0) return;
  for(int i=0;i<labels->numBins("");i++){
    labels->setIndex(i);
    lblVec->push_back(labels->getLabel());
    std::cout << lblVec->back() <<std::endl;
  }
}

void MakeSpinFits::MakeSignalFitForFit(TString tag, TString mcName){
  std::cout << "\n\n\nMakeSpinFits::MakeSignalFitForFit\n\n\n" << std::endl;
  if(ws==0) return;

  RooRealVar mass = *(ws->var("mass"));

  //instantiate parameters copy from original signal fits
  RooRealVar mean( Form("Data_%s_FIT_mean",mcName.Data()), "", meanStart,meanStart-0.5,meanStart+0.5);
  RooRealVar sig1( Form("Data_%s_FIT_%s_sigma1",mcName.Data(),tag.Data()), "", 
		   ws->var( Form("%s_FIT_%s_sigma1",mcName.Data(),tag.Data()))->getVal());
  RooRealVar sig2( Form("Data_%s_FIT_%s_sigma2",mcName.Data(),tag.Data()), "", 
		   ws->var( Form("%s_FIT_%s_sigma2",mcName.Data(),tag.Data()))->getVal());
  RooRealVar sig3( Form("Data_%s_FIT_%s_sigma3",mcName.Data(),tag.Data()), "", 
		   ws->var( Form("%s_FIT_%s_sigma3",mcName.Data(),tag.Data()))->getVal());
  RooRealVar f1( Form("Data_%s_FIT_%s_f1",mcName.Data(),tag.Data()), "", 
		   ws->var( Form("%s_FIT_%s_f1",mcName.Data(),tag.Data()))->getVal());
  RooRealVar f2( Form("Data_%s_FIT_%s_f2",mcName.Data(),tag.Data()), "", 
		   ws->var( Form("%s_FIT_%s_f2",mcName.Data(),tag.Data()))->getVal());

  RooRealVar *sigCB=0, *alphaCB=0, *nCB=0, *fCB=0;
  if(useCB){
    sigCB = new RooRealVar(Form("Data_%s_FIT_%s_sigmaCB",mcName.Data(),tag.Data()),"",
			   ws->var(Form("%s_FIT_%s_sigmaCB",mcName.Data(),tag.Data()))->getVal());

    alphaCB = new RooRealVar(Form("Data_%s_FIT_%s_alphaCB",mcName.Data(),tag.Data()),"",
			     ws->var(Form("%s_FIT_%s_alphaCB",mcName.Data(),tag.Data()))->getVal());
    nCB = new RooRealVar(Form("Data_%s_FIT_%s_nCB",mcName.Data(),tag.Data()),"",
			 ws->var(Form("%s_FIT_%s_nCB",mcName.Data(),tag.Data()))->getVal());
    fCB = new RooRealVar(Form("Data_%s_FIT_%s_fCB",mcName.Data(),tag.Data()),"",
		     ws->var(Form("%s_FIT_%s_fCB",mcName.Data(),tag.Data()))->getVal());
    sigCB->setConstant(kTRUE);
    alphaCB->setConstant(kTRUE);
    nCB->setConstant(kTRUE);
    fCB->setConstant(kTRUE);
  }
  //mean.setConstant(kTRUE);
  sig1.setConstant(kTRUE);
  sig2.setConstant(kTRUE);
  sig2.setConstant(kTRUE);
  f1.setConstant(kTRUE);
  f2.setConstant(kTRUE);


  //triple gaussian
 RooGaussian g1( Form("Data_%s_FIT_%s_g1",mcName.Data(),tag.Data()), "",mass,mean,sig1); 
  RooGaussian g2( Form("Data_%s_FIT_%s_g2",mcName.Data(),tag.Data()), "",mass,mean,sig2); 
  RooGaussian g3( Form("Data_%s_FIT_%s_g3",mcName.Data(),tag.Data()), "",mass,mean,sig3); 
  
  RooCBShape *cb=0;
  if(useCB) cb = new RooCBShape(Form("Data_%s_FIT_%s_cb",mcName.Data(),tag.Data()),"",mass,mean,*sigCB,*alphaCB,*nCB);

  RooAddPdf *SignalModel=0;
  if(useCB){
    SignalModel = new RooAddPdf(Form("Data_%s_FIT_%s",mcName.Data(),tag.Data()),"Signal Model",RooArgList(g1,g2,g3,*cb),RooArgList(f1,f2,*fCB));
  }else{
    SignalModel = new RooAddPdf(Form("Data_%s_FIT_%s",mcName.Data(),tag.Data()),"Signal Model",RooArgList(g1,g2,g3),RooArgList(f1,f2));
  }

  //RooAddPdf SignalModel( Form("Data_%s_FIT_%s",mcName.Data(),tag.Data()),"Signal Model",RooArgList(g1,g2,g3),RooArgList(f1,f2));
  
  ws->import(*SignalModel);
  delete SignalModel;
  if(useCB){
    delete sigCB;
    delete alphaCB;
    delete nCB;
    delete fCB;
    delete cb;
  }
}

void MakeSpinFits::MakeSignalFit(TString tag, TString mcName,float cosTlow, float cosThigh){
  RooTrace::active(__DO_TRACE);
  //RooTrace::verbose(kTRUE);
  std::cout << "\n\n\nMakeSpinFits::MakeSignalFit\n\n\n" << ws << std::endl;
  if(ws==0) return;
  TString inputTag = Form("%s_%s",mcName.Data(),tag.Data());
  TString outputTag = Form("%s_FIT_%s",mcName.Data(),tag.Data());
  ws->Print();
  RooRealVar mass = *(ws->var("mass"));
  RooRealVar cosT = *(ws->var("cosT"));
  cosT.setBins(10);

  if(cosTlow > -1 || cosThigh < 1){
    outputTag+=Form("_cosT_%0.2f_%0.2f",cosTlow,cosThigh);
  }

  //signal fit parameters -- Triple Gaussian
  RooRealVar mean(Form("%s_mean",outputTag.Data()),Form("%s_mean",outputTag.Data()),meanStart,meanLow,meanHigh);
  RooRealVar sig1(Form("%s_sigma1",outputTag.Data()),Form("%s_sigma1",outputTag.Data()),1,0.1,6);
  RooRealVar sig2(Form("%s_sigma2",outputTag.Data()),Form("%s_sigma2",outputTag.Data()),4,2,10);
  RooRealVar sig3(Form("%s_sigma3",outputTag.Data()),Form("%s_sigma3",outputTag.Data()),7,5,40);
  RooRealVar f1(Form("%s_f1",outputTag.Data()),Form("%s_f1",outputTag.Data()),0.1,0.00011,1);
  RooRealVar f2(Form("%s_f2",outputTag.Data()),Form("%s_f2",outputTag.Data()),0.1,0.0001,1);
  RooGaussian g1(Form("%s_g1",outputTag.Data()),Form("%s_g1",outputTag.Data()),mass,mean,sig1);
  RooGaussian g2(Form("%s_g2",outputTag.Data()),Form("%s_g2",outputTag.Data()),mass,mean,sig2);
  RooGaussian g3(Form("%s_g3",outputTag.Data()),Form("%s_g3",outputTag.Data()),mass,mean,sig3);
  //RooAddPdf SignalModel(outputTag.Data(),"Signal Model",RooArgList(g1,g2,g3),RooArgList(f1,f2));

  //signal fit parameters -- CB times triple gaussian
  RooRealVar sigCB(Form("%s_sigmaCB",outputTag.Data()),Form("%s_sigmaCB",outputTag.Data()),1,0.2,10);
  RooRealVar alphaCB(Form("%s_alphaCB",outputTag.Data()),Form("%s_alphaCB",outputTag.Data()),1,0.,1000);
  RooRealVar nCB(Form("%s_nCB",outputTag.Data()),Form("%s_nCB",outputTag.Data()),1,0.,1000);
  RooRealVar fCB(Form("%s_fCB",outputTag.Data()),Form("%s_fCB",outputTag.Data()),0.1,0,1);
  
  RooCBShape cb(Form("%s_cb",outputTag.Data()),Form("%s_cb",outputTag.Data()),mass,mean,sigCB,alphaCB,nCB);

  RooAddPdf *SignalModel=0;
  if(useCB){
    SignalModel = new RooAddPdf(outputTag.Data(),"Signal Model",RooArgList(g1,g2,g3,cb),RooArgList(f1,f2,fCB));
  }else{
    SignalModel = new RooAddPdf(outputTag.Data(),"Signal Model",RooArgList(g1,g2,g3),RooArgList(f1,f2));
  }

   
  
  //signal MC data
  RooDataSet *ds=0;
  if(tag!="Combined") ds = (RooDataSet*)ws->data( Form("%s_Combined",mcName.Data()) )->reduce(TString("evtcat==evtcat::")+tag.Data());
  else ds = (RooDataSet*)ws->data( Form("%s_Combined",mcName.Data()) );

  if(cosTlow > -1 || cosThigh < 1){
    ds = (RooDataSet*)ds->reduce( Form("cosT >= %0.3f && cosT < %0.3f",cosTlow,cosThigh) );
  }

  RooDataHist hist(Form("%s_cosThist",outputTag.Data()),"Data Hist for cos(theta)",RooArgSet(cosT),*ds);
  RooHistPdf cosTkde(Form("%s_cosTpdf",outputTag.Data()),"Hist PDF for cos(theta)",RooArgSet(cosT),hist);

  SignalModel->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU));
  RooFitResult *res = SignalModel->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU));
  //SignalModelCB.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU));
  //RooFitResult *resCB = SignalModelCB.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU));
  std::cout << res <<std::endl;
  res->SetName(Form("%s_FitResult",outputTag.Data()));
  //resCB->SetName(Form("%s_FitResultCB",outputTag.Data()));

  cosTkde.fitTo(*ds);

  //set all parameters as constant
  mean.setConstant(kTRUE);
  sig1.setConstant(kTRUE);
  sig2.setConstant(kTRUE);
  sig3.setConstant(kTRUE);
  f1.setConstant(kTRUE);
  f2.setConstant(kTRUE);

  //compute the full width at half maximum and the sigma effective
  float FWHM = computeFWHM(SignalModel,mean.getVal(),&mass);
  float sige = computeSigEff(SignalModel,mean.getVal(),&mass);
  RooRealVar rvFWHM(Form("%s_FWHM",outputTag.Data()),"",FWHM);


  RooRealVar rvSE(Form("%s_sigmaEff",outputTag.Data()),"",sige);
  
  ws->import(*SignalModel);
  ws->import(cosTkde);
  ws->import(rvFWHM);
  ws->import(rvSE);
  ws->import(*res);

  delete SignalModel;
  //delete ds;
}

float MakeSpinFits::computeFWHM(RooAbsPdf* pdf, float mean,RooRealVar *var){
  std::cout << "MakeSpinFits::computeFWHM" << std::endl;
  //computes the FWHM for a PDF
  //ONLY VALID FOR SYMMETRIC PDFS!!!
  //need to fix to make more general
  var->setVal(mean);
  float peakV = pdf->getVal();
  
  const float step = 0.05;

  float firstNegOffset=0;
  for(float offset = 0; offset < step*200; offset+=step){
    var->setVal(mean+offset);
    float thisV = pdf->getVal();
    float ratioDiff = thisV/peakV-0.5; //wait for this to go negative
    std::cout << offset << "    " << ratioDiff;
    if(ratioDiff <= 0){
      firstNegOffset = offset;
      break;
    }
  }

  float diffNeg = fabs( 0.5 - pdf->getVal()/peakV);

  var->setVal(firstNegOffset-step);
  float diffPos = fabs( 0.5 - pdf->getVal()/peakV);  

  if(diffPos<diffNeg) return 2*firstNegOffset - 2*step;
  return 2*firstNegOffset;
}

float MakeSpinFits::computeSigEff(RooAbsPdf* pdf, float mean, RooRealVar *var){
  std::cout << "MakeSpinFits::computeSigEff" << std::endl;
  // compute the sigma_eff for the pdf
  //only valid for symmetric pdfs

  float width=0.0;

  float below=0, belowCov=0;
  float above=0, aboveCov=0;

  while(true){
    width+=0.01;
    var->setRange("sigEff",mean-width,mean+width);
    RooAbsReal * integral = pdf->createIntegral(*var,RooFit::NormSet(*var),RooFit::Range("sigEff"));
    float cov = integral->getVal();
    delete integral;
    if(cov > 0.683){
      above=width;
      aboveCov=cov;
      break;
    }else{
      below=width;
      belowCov=cov;
    }
  }

  return (above*(aboveCov-0.683)+below*(0.683-belowCov))/(aboveCov-belowCov); //weighted average

}

void MakeSpinFits::MakeCombinedSignalSpin(TString mcName){
  std::cout << "MakeSpinFits::MakeCombinedSignalSpin" <<std::endl;
  RooRealVar cosT = *(ws->var("cosT"));

  RooDataSet *data = (RooDataSet*)ws->data( Form("%s_Combined",mcName.Data()) );

  RooDataHist hist(Form("%s_FIT_cosThist", mcName.Data()),"",cosT,*data);
  RooHistPdf  pdf (Form("%s_FIT_cosTpdf",  mcName.Data()),"",cosT, hist);
  ws->import(pdf);
}

void MakeSpinFits::AddCombinedBkgOnlySWeight(TString mcName){
  std::cout << "MakeSpinFits::AddCombinedBkgOnlySWeight" <<std::endl;
  if(ws==0) return;

  RooRealVar *mass = (RooRealVar*)ws->var("mass");

  std::vector<TString>::const_iterator catIt = catLabels.begin();
  for(; catIt != catLabels.end(); catIt++){
    cout << *catIt <<endl;

    RooDataSet * data = (RooDataSet*)ws->data("Data_Combined")->reduce(TString("evtcat==evtcat::")+*catIt);
    RooAbsPdf *signalModel = ws->pdf( Form("%s_FIT_%s",mcName.Data(),catIt->Data()) );
    RooAbsPdf *bkgModel    = ws->pdf( Form("Data_BKGFIT_%s_bkgModel",catIt->Data() ) );    

    RooFormulaVar* Nsig = (RooFormulaVar*)ws->obj( Form("Data_%s_FULLFIT_%s_Nsig",mcName.Data(),catIt->Data()) );
    RooFormulaVar* Nbkg = (RooFormulaVar*)ws->obj( Form("Data_%s_FULLFIT_%s_Nbkg",mcName.Data(),catIt->Data()) );

    cout << data << "   " << signalModel << "  " << bkgModel << "  " << Nsig << "  " << Nbkg <<endl;

    //instantiate the SPlot class
    MakeSpinSPlot splotter(data);
    splotter.addSpecies("background",bkgModel,Nbkg->getVal());
    splotter.addSpecies("signal",signalModel,Nsig->getVal());

    splotter.addVariable(mass);
    cout << "Calculating SWeight ... " <<endl;
    splotter.calculate();
    cout << "Done" <<endl;

    RooDataSet * sweights = splotter.getSWeightDataSet();
    sweights->merge(data);
    sweights->SetName( Form("Data_%s_%s_WithSWeight",catIt->Data(),mcName.Data()) );
    
    //get list of variables in the original dataset
    RooArgSet dsVars;
    const RooArgSet *orig = data->get(0);
    RooLinkedListIter it = orig->iterator();
    while(it.Next()){
      std::cout << (*it)->GetName() << std::endl;
      if(ws->var( (*it)->GetName()) ==0 ) continue;
      dsVars.add(*(ws->var( (*it)->GetName()) ));
    }

    RooArgSet *swVars = splotter.getSWeightVars();
    dsVars.add(*(RooRealVar*)(swVars->find("signal_sw")) );
    dsVars.add(*(RooRealVar*)(swVars->find("background_sw")) );
    //make weighted datasets
    RooDataSet sig(Form("Data_%s_%s_sigWeight",catIt->Data(),mcName.Data()),"",sweights,dsVars,0,"signal_sw");
    RooDataSet bkg(Form("Data_%s_%s_bkgWeight",catIt->Data(),mcName.Data()),"",sweights,dsVars,0,"background_sw");
    

    ws->import(*sweights);
    ws->import(sig);
    ws->import(bkg);
    delete data;
  }    
}

void MakeSpinFits::MakeFullSBFit(TString mcName,bool cosTBinned){
  std::cout << "MakeSpinFits::MakeFullSBFit" <<std::endl;
  if(ws==0) return;

  RooRealVar *mass = ws->var("mass");

  RooCategory *cat;
  if(cosTBinned) cat = ((RooCategory*)ws->obj("evtcat_cosT"));
  else cat = ((RooCategory*)ws->obj("evtcat"));

  TString fitTag = "FULLSBFIT";
  if(cosTBinned) fitTag  = "FULLSBCOSTFIT";

  RooSimultaneous *combFit = new RooSimultaneous(Form("Data_%s_%s",mcName.Data(),fitTag.Data()),"",*cat);

  RooRealVar *nSig = new RooRealVar( Form("Data_%s_%s_Nsig",mcName.Data(),fitTag.Data() ), "", 300, 0,1e6);
  
  std::vector<TString> *cats = (cosTBinned ? &cosTcatLabels : &catLabels);

  std::vector<TString>::const_iterator catIt = cats->begin();
  for(; catIt != cats->end(); catIt++){
    cout << *catIt <<endl;
    TString dsName = "_Combined";
    if(cosTBinned) dsName = "_Combined_CosTBin";
    TString selectionString = TString("evtcat==evtcat::")+*catIt;
    if(cosTBinned) selectionString = TString("evtcat_cosT==evtcat_cosT::")+*catIt;

    //get the signal Model
    RooAbsPdf *signalModel = ws->pdf( Form("Data_%s_FIT_%s",mcName.Data(),catIt->Data()) );
    //background model
    RooAbsPdf* BkgShape;
    
    switch(fitType){
    case kExp:
      {
	//double exponential
	RooRealVar* alpha1 = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_alpha1",mcName.Data(),catIt->Data()),"alpha1",-0.1,-1.,0.);
	RooRealVar* alpha2 = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_alpha2",mcName.Data(),catIt->Data()),"alpha2",-0.1,-1.,0.);
	RooRealVar* f_bkg  = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_f",mcName.Data(),catIt->Data()),"f_bkg",0.1,0,1);
	RooExponential* exp1 = new RooExponential(Form("Data_%s_FULLSBFIT_%s_exp1",mcName.Data(),catIt->Data()),"exp1",*mass,*alpha1);
	RooExponential* exp2 = new RooExponential(Form("Data_%s_FULLSBFIT_%s_exp2",mcName.Data(),catIt->Data()),"exp2",*mass,*alpha2);
	
	BkgShape = new RooAddPdf(Form("Data_%s_FULLSBFIT_%s_bkgShape",mcName.Data(),catIt->Data()),"Background Model",
				 RooArgList(*exp1,*exp2),*f_bkg);
	break;
      }
    case kPoly:
      {
	//5th order polynomial
	RooRealVar *pC = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_pC",mcName.Data(),catIt->Data()),"pC",1);
	RooRealVar *p0 = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_p0",mcName.Data(),catIt->Data()),"p0",0,-10,10);
	RooRealVar *p1 = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_p1",mcName.Data(),catIt->Data()),"p1",0,-10,10);
	RooRealVar *p2 = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_p2",mcName.Data(),catIt->Data()),"p2",0,-10,10);
	RooRealVar *p3 = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_p3",mcName.Data(),catIt->Data()),"p3",0,-10,10);
	RooRealVar *p4 = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_p4",mcName.Data(),catIt->Data()),"p4",0,-10,10);
	//enforce all coefficients positive
	RooFormulaVar *pCmod = new RooFormulaVar(Form("Data_%s_FULLSBFIT_%s_pCmod",mcName.Data(),catIt->Data()),"","@0*@0",*pC);
	RooFormulaVar *p0mod = new RooFormulaVar(Form("Data_%s_FULLSBFIT_%s_p0mod",mcName.Data(),catIt->Data()),"","@0*@0",*p0);
	RooFormulaVar *p1mod = new RooFormulaVar(Form("Data_%s_FULLSBFIT_%s_p1mod",mcName.Data(),catIt->Data()),"","@0*@0",*p1);
	RooFormulaVar *p2mod = new RooFormulaVar(Form("Data_%s_FULLSBFIT_%s_p2mod",mcName.Data(),catIt->Data()),"","@0*@0",*p2);
	RooFormulaVar *p3mod = new RooFormulaVar(Form("Data_%s_FULLSBFIT_%s_p3mod",mcName.Data(),catIt->Data()),"","@0*@0",*p3);
	RooFormulaVar *p4mod = new RooFormulaVar(Form("Data_%s_FULLSBFIT_%s_p4mod",mcName.Data(),catIt->Data()),"","@0*@0",*p4);
	
  

	RooArgList *args;
	if(cosTBinned){
	  args = new RooArgList(*pCmod,*p0mod,*p1mod,*p2mod,*p3mod);
	}else{
	  args = new RooArgList(*pCmod,*p0mod,*p1mod,*p2mod,*p3mod,*p4mod);
	}
	
	BkgShape = new RooBernstein(Form("Data_%s_FULLSBFIT_%s_bkgShape",mcName.Data(),catIt->Data()),"Background Model",*mass,*args);
	break;
      }
    default:
      std::cout << "INVALID BACKGROUND MODEL" << std::endl;
      assert(false);
      break;
    }

    double thisCatB     = ws->data("Data"+dsName)->sumEntries(selectionString);
    RooRealVar *Nbkg = new RooRealVar(Form("Data_%s_FULLSBFIT_%s_Nbkg",mcName.Data(),catIt->Data()),"N background Events",thisCatB,0,1e9);
    //extended fit model
    RooExtendPdf *bkgModel = new RooExtendPdf(Form("Data_%s_FULLSBFIT_%s_bkgModel",mcName.Data(),catIt->Data()),"Background Model",*BkgShape,*Nbkg);


    //compute the fraction of events expected in this category
    double totalEvents  = ws->data(mcName+dsName)->sumEntries(); 
    double thisCat =  ws->data(mcName+dsName)->sumEntries(selectionString);
    double thisFrac = thisCat/totalEvents;
    double thisFracE = thisFrac * TMath::Sqrt(1/thisCat+1/totalEvents);

    RooRealVar *fSig    = new RooRealVar( Form("Data_%s_%s_%s_fsig",mcName.Data(),fitTag.Data(),catIt->Data() ), "", thisFrac);
    //fix the signal fraction
    RooFormulaVar *thisNsig = new RooFormulaVar(Form("Data_%s_%s_%s_Nsig",mcName.Data(),fitTag.Data(),catIt->Data() ),"","@0*@1",RooArgSet(*nSig,*fSig));

    //RooExtendPdf *exSignalModel = new RooExtendPdf(Form("Data_%s_FULLSBFIT_%s_signalModel",mcName.Data(),catIt->Data()),"",*signalModel,*thisNsig);

    //build the combined fit model
    RooAddPdf *comb = new RooAddPdf(Form("Data_%s_%s_%s",mcName.Data(),fitTag.Data(),catIt->Data()),"",RooArgList(*signalModel,*bkgModel),
				    RooArgList(*thisNsig,*Nbkg) );

    combFit->addPdf(*comb,*catIt);
  }    

  RooDataSet *ds=0;
  if(cosTBinned) ds = (RooDataSet*)ws->data("Data_Combined_CosTBin");
  else ds = (RooDataSet*)ws->data("Data_Combined");

  combFit->fitTo(*ds,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=combFit->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_%s_%s_fitResult",mcName.Data(),fitTag.Data()) );
  
  ws->import(*combFit);
  ws->import(*res);

}


void MakeSpinFits::MakeCombinedSignalTest(TString mcName,bool cosTBinned){
  std::cout << "MakeSpinFits::MakeCombinedSignalTest" <<std::endl;
  if(ws==0) return;

  RooCategory *cat;
  if(cosTBinned) cat = ((RooCategory*)ws->obj("evtcat_cosT"));
  else cat = ((RooCategory*)ws->obj("evtcat"));

  TString fitTag = "FULLFIT";
  if(cosTBinned) fitTag  = "FULLCOSTFIT";

  RooSimultaneous *combFit = new RooSimultaneous(Form("Data_%s_%s",mcName.Data(),fitTag.Data()),"",*cat);

  RooRealVar *nSig = new RooRealVar( Form("Data_%s_%s_Nsig",mcName.Data(),fitTag.Data() ), "", 300, 0,1e6);

  std::vector<TString> *cats = (cosTBinned ? &cosTcatLabels : &catLabels);

  std::vector<TString>::const_iterator catIt = cats->begin();
  for(; catIt != cats->end(); catIt++){
    cout << *catIt <<endl;
    TString dsName = "_Combined";
    if(cosTBinned) dsName = "_Combined_CosTBin";
    TString selectionString = TString("evtcat==evtcat::")+*catIt;
    if(cosTBinned) selectionString = TString("evtcat_cosT==evtcat_cosT::")+*catIt;

    //get the signal Model
    RooAbsPdf *signalModel = ws->pdf( Form("Data_%s_FIT_%s",mcName.Data(),catIt->Data()) );
    RooAbsPdf *bkgModel    = ws->pdf( Form("Data_BKGFIT_%s_bkgModel",catIt->Data() ) );
    std::cout << signalModel << "   " << bkgModel << std::endl;

    double thisCatB     = ws->data("Data"+dsName)->sumEntries(selectionString);

    RooRealVar *nBkg = new RooRealVar( Form("Data_%s_%s_%s_Nbkg",mcName.Data(),fitTag.Data(),catIt->Data() ), "", thisCatB, 0,1e9);
    //    RooRealVar *fBkg = new RooRealVar( Form("Data_%s_%s_%s_fbkg",mcName.Data(),fitTag.Data(),catIt->Data()), "",
    //thisCatB/totalEventsB,thisCatB/totalEventsB-2*TMath::Sqrt(thisCatB)/totalEventsB,thisCatB/totalEventsB+2*TMath::Sqrt(thisCatB)/totalEventsB); 
    //RooFormulaVar *thisNbkg = new RooFormulaVar(Form("Data_%s_%s_%s_Nbkg",mcName.Data(),fitTag.Data(),catIt->Data() ),"","@0*@1",RooArgSet(*nBkg,*fBkg));

    //compute the fraction of events expected in this category
    double totalEvents  = ws->data(mcName+dsName)->sumEntries(); 
    double thisCat =  ws->data(mcName+dsName)->sumEntries(selectionString);
    double thisFrac = thisCat/totalEvents;
    double thisFracE = thisFrac * TMath::Sqrt(1/thisCat+1/totalEvents);

    RooRealVar *fSig    = new RooRealVar( Form("Data_%s_%s_%s_fsig",mcName.Data(),fitTag.Data(),catIt->Data() ), "", thisFrac);
    //fix the signal fraction
    RooFormulaVar *thisNsig = new RooFormulaVar(Form("Data_%s_%s_%s_Nsig",mcName.Data(),fitTag.Data(),catIt->Data() ),"","@0*@1",RooArgSet(*nSig,*fSig));

    //RooExtendPdf *exSignalModel = new RooExtendPdf(Form("Data_%s_FULLFIT_%s_signalModel",mcName.Data(),catIt->Data()),"",*signalModel,*thisNsig);

    //build the combined fit model
    RooAddPdf *comb = new RooAddPdf(Form("Data_%s_%s_%s",mcName.Data(),fitTag.Data(),catIt->Data()),"",RooArgList(*signalModel,*bkgModel),
				    RooArgList(*thisNsig,*nBkg) );

    combFit->addPdf(*comb,*catIt);
  }    

  RooDataSet *ds=0;
  if(cosTBinned) ds = (RooDataSet*)ws->data("Data_Combined_CosTBin");
  else ds = (RooDataSet*)ws->data("Data_Combined");

  combFit->fitTo(*ds,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=combFit->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_%s_%s_fitResult",mcName.Data(),fitTag.Data()) );
  
  ws->import(*nSig);
  ws->import(*combFit);
  ws->import(*res);

}

void MakeSpinFits::MakeFloatingSignalTest(TString mcName){
  std::cout << "MakeSpinFits::MakeFloatingSignalTest" <<std::endl;
  if(ws==0) return;

  RooCategory *cat = ((RooCategory*)ws->obj("evtcat"));

  RooSimultaneous *combFit = new RooSimultaneous(Form("Data_%s_INDFIT",mcName.Data()),"",*cat);

  std::vector<TString>::const_iterator catIt = catLabels.begin();
  for(; catIt != catLabels.end(); catIt++){
    cout << *catIt <<endl;
    
    RooAbsPdf *signalModel = ws->pdf( Form("Data_%s_FIT_%s",mcName.Data(),catIt->Data()) );
    RooAbsPdf *bkgModel    = ws->pdf( Form("Data_BKGFIT_%s_bkgModel",catIt->Data() ) );
    std::cout << signalModel << "   " << bkgModel << std::endl;

    double totalEventsB = ws->data("Data_Combined")->sumEntries();
    double thisCatB     = ws->data("Data_Combined")->sumEntries(TString("evtcat==evtcat::")+*catIt);

    RooRealVar *nBkg = new RooRealVar( Form("Data_%s_INDFIT_%s_Nbkg",mcName.Data(),catIt->Data()), "", thisCatB,0,1e9 ); 

    //signal yield floated in each category
    RooRealVar *nSig    = new RooRealVar( Form("Data_%s_INDFIT_%s_Nsig",mcName.Data(),catIt->Data() ), "", 0,0,1e9);

    //RooExtendPdf *exSignalModel = new RooExtendPdf(Form("Data_%s_INDFIT_%s_signalModel",mcName.Data(),catIt->Data()),"",*signalModel,*thisNsig);

    RooAddPdf *comb = new RooAddPdf(Form("Data_%s_INDFIT_%s",mcName.Data(),catIt->Data()),"",RooArgList(*signalModel,*bkgModel),
				    RooArgList(*nSig,*nBkg) );

    combFit->addPdf(*comb,*catIt);
  }    

  RooDataSet *ds = (RooDataSet*)ws->data("Data_Combined");

  combFit->fitTo(*ds,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=combFit->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_%s_INDFIT_fitResult",mcName.Data()) );
  
  ws->import(*combFit);
  ws->import(*res);

}

void MakeSpinFits::Make2DCombinedSignalTest(TString massMcName,TString costMcName){
  std::cout << "MakeSpinFits::Make2DCombinedSignalTest" <<std::endl;
  if(ws==0) return;

  RooCategory *cat = ((RooCategory*)ws->obj("evtcat"));

  RooSimultaneous *combFit = new RooSimultaneous(Form("Data_m_%s_c_%s_FULL2DFIT",massMcName.Data(),costMcName.Data()),"",*cat);

  RooRealVar *nSig = new RooRealVar( Form("Data_m_%s_c_%s_FULL2DFIT_Nsig",massMcName.Data(),costMcName.Data() ), "", 0, 0,1e9);
  RooRealVar *nBkg = new RooRealVar( Form("Data_m_%s_c_%s_FULL2DFIT_Nbkg",massMcName.Data(),costMcName.Data() ), "", 0, 0,1e9);

  RooRealVar *mass = ws->var("mass");
  RooRealVar *cosT = ws->var("cosT");

  std::vector<TString>::const_iterator catIt = catLabels.begin();
  for(; catIt != catLabels.end(); catIt++){
    cout << *catIt <<endl;
    
    //mass models
    RooAbsPdf *signalMassModel = ws->pdf( Form("Data_%s_FIT_%s",massMcName.Data(),catIt->Data()) );
    RooAbsPdf *bkgMassModel    = ws->pdf( Form("Data_BKGFIT_%s_bkgShape",catIt->Data() ) );

    //cos(theta) models
    RooAbsPdf *signalCosModel = ws->pdf(Form("%s_FIT_%s_cosTpdf",costMcName.Data(),catIt->Data()));
    RooAbsPdf *bkgCosModel    = ws->pdf(Form("Data_BKGFIT_%s_cosTpdf",catIt->Data()));

    std::cout << bkgCosModel <<std::endl;

    double totalEventsB = ws->data("Data_Combined")->sumEntries();
    double thisCatB     = ws->data("Data_Combined")->sumEntries(TString("evtcat==evtcat::")+*catIt);

    RooRealVar *fBkg = new RooRealVar( Form("Data_m_%s_c_%s_FULL2DFIT_%s_fbkg",massMcName.Data(),costMcName.Data(),catIt->Data()), "", thisCatB/totalEventsB ); 
    RooFormulaVar *thisNbkg = new RooFormulaVar(Form("Data_m_%s_c_%s_FULL2DFIT_%s_Nbkg",massMcName.Data(),costMcName.Data(),catIt->Data() ),"","@0*@1",RooArgSet(*nBkg,*fBkg));


    double totalEvents  = ws->data(massMcName+"_Combined")->sumEntries();
    double thisCat =  ws->data(massMcName+"_Combined")->sumEntries(TString("evtcat==evtcat::")+*catIt);
    double thisFrac = thisCat/totalEvents;
    double thisFracE = thisFrac * TMath::Sqrt(1/thisCat+1/totalEvents);
    //fix signal fraction
    RooRealVar *fSig    = new RooRealVar( Form("Data_m_%s_c_%s_FULL2DFIT_%s_fsig",massMcName.Data(),costMcName.Data(),catIt->Data() ), "", thisFrac);
    RooFormulaVar *thisNsig = new RooFormulaVar(Form("Data_m_%s_c_%s_FULL2DFIT_%s_Nsig",massMcName.Data(),costMcName.Data(),catIt->Data() ),"","@0*@1",RooArgSet(*nSig,*fSig));

    //RooExtendPdf *exSignalModel = new RooExtendPdf(Form("Data_%s_FULL2DFIT_%s_signalModel",mcName.Data(),catIt->Data()),"",*signalModel,*thisNsig);

    /*
    std::cout << "sig model MASS integral: " << signalMassModel->createIntegral(*mass)->getVal() <<std::endl;
    std::cout << "bkg model MASS integral: " << bkgMassModel->createIntegral(*mass)->getVal() <<std::endl;

    std::cout << "sig model COST integral: " << signalCosModel->createIntegral(*cosT)->getVal() <<std::endl;
    std::cout << "bkg model COST integral: " << bkgCosModel->createIntegral(*cosT)->getVal() <<std::endl;
    */
    //build 2D signal/bkg models
    RooProdPdf *signalModel = new RooProdPdf(Form("Data_m_%s_c_%s_FULL2DFIT_%s_signalModel",massMcName.Data(),costMcName.Data(),catIt->Data()),"",
					     RooArgList(*signalMassModel,*signalCosModel));
    RooProdPdf *bkgModel = new RooProdPdf(Form("Data_m_%s_c_%s_FULL2DFIT_%s_bkgModel",massMcName.Data(),costMcName.Data(),catIt->Data()),"",
					     RooArgList(*bkgMassModel,*bkgCosModel));
    /*
    std::cout << "sig model integral: " << signalModel->createIntegral(RooArgList(*mass,*cosT))->getVal() <<std::endl;
    std::cout << "bkg model integral: " << bkgModel->createIntegral(RooArgList(*mass,*cosT))->getVal() <<std::endl;
    */
    //build combined fit model
    RooAddPdf *comb = new RooAddPdf(Form("Data_m_%s_c_%s_FULL2DFIT_%s",massMcName.Data(),costMcName.Data(),catIt->Data()),"",RooArgList(*signalModel,*bkgModel),
				    RooArgList(*thisNsig,*thisNbkg) );

    combFit->addPdf(*comb,*catIt);
  }    

  RooDataSet *ds = (RooDataSet*)ws->data("Data_Combined");

  combFit->fitTo(*ds,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=combFit->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_m_%s_c_%s_FULL2DFIT_fitResult",massMcName.Data(),costMcName.Data()) );
  
  ws->import(*combFit);
  ws->import(*res);

}

void MakeSpinFits::Make2DTemplateSignalTest(TString mcName){
  std::cout << "MakeSpinFits::Make2DCombinedSignalTest" <<std::endl;
  if(ws==0) return;

  RooCategory *cat = ((RooCategory*)ws->obj("evtcat"));

  RooSimultaneous *combFit = new RooSimultaneous(Form("Data_%s_TEMPLATE2DFIT",mcName.Data()),"",*cat);

  RooRealVar *nSig = new RooRealVar( Form("Data_%s_TEMPLATE2DFIT_Nsig",mcName.Data() ), "", 0, 0,1e9);
  RooRealVar *nBkg = new RooRealVar( Form("Data_%s_TEMPLATE2DFIT_Nbkg",mcName.Data() ), "", 0, 0,1e9);

  RooRealVar *mass = ws->var("mass");
  RooRealVar *cosT = ws->var("cosT");

  std::vector<RooAbsPdf*>  pdfsToDelete;

  std::vector<TString>::const_iterator catIt = catLabels.begin();
  for(; catIt != catLabels.end(); catIt++){
    cout << *catIt <<endl;
    
    //mass models
    RooAbsData *signalData = ws->data(Form("%s_Combined",mcName.Data()))->reduce(TString("evtcat==evtcat::"+*catIt));
    RooDataHist *signalDataHist = new RooDataHist( Form("Data_%s_TEMPLATE2DFIT_%s_signalDataHist",mcName.Data(),catIt->Data()), "", RooArgSet(*mass,*cosT),*signalData);
    RooHistPdf *signalHistPdf = new RooHistPdf( Form("Data_%s_TEMPLATE2DFIT_%s_signalHistPdf",mcName.Data(),catIt->Data()), "",  RooArgSet(*mass,*cosT), * signalDataHist);

    pdfsToDelete.push_back(signalHistPdf);

    RooAbsPdf *bkgMassModel    = ws->pdf( Form("Data_BKGFIT_%s_bkgShape",catIt->Data() ) );

    //cos(theta) models
    RooAbsPdf *bkgCosModel    = ws->pdf(Form("Data_BKGFIT_%s_cosTpdf",catIt->Data()));

    double totalEventsB = ws->data("Data_Combined")->sumEntries();
    double thisCatB     = ws->data("Data_Combined")->sumEntries(TString("evtcat==evtcat::")+*catIt);

    RooRealVar *fBkg = new RooRealVar( Form("Data_%s_TEMPLATE2DFIT_%s_fbkg",mcName.Data(),catIt->Data()), "", thisCatB/totalEventsB ); 
    RooFormulaVar *thisNbkg = new RooFormulaVar(Form("Data_%s_TEMPLATE2DFIT_%s_Nbkg",mcName.Data(),catIt->Data() ),"","@0*@1",RooArgSet(*nBkg,*fBkg));


    double totalEvents  = ws->data(mcName+"_Combined")->sumEntries();
    double thisCat =  ws->data(mcName+"_Combined")->sumEntries(TString("evtcat==evtcat::")+*catIt);
    double thisFrac = thisCat/totalEvents;
    double thisFracE = thisFrac * TMath::Sqrt(1/thisCat+1/totalEvents);
    //fix signal fraction
    RooRealVar *fSig    = new RooRealVar( Form("Data_%s_TEMPLATE2DFIT_%s_fsig",mcName.Data(),catIt->Data() ), "", thisFrac);
    RooFormulaVar *thisNsig = new RooFormulaVar(Form("Data_%s_TEMPLATE2DFIT_%s_Nsig",mcName.Data(),catIt->Data() ),"","@0*@1",RooArgSet(*nSig,*fSig));

    //RooExtendPdf *exSignalModel = new RooExtendPdf(Form("Data_%s_TEMPLATE2DFIT_%s_signalModel",mcName.Data(),catIt->Data()),"",*signalModel,*thisNsig);


    //build 2D signal/bkg models
    RooProdPdf *bkgModel = new RooProdPdf(Form("Data_%s_TEMPLATE2DFIT_%s_bkgModel",mcName.Data(),catIt->Data()),"",
					     RooArgList(*bkgMassModel,*bkgCosModel));
    
    //std::cout << "bkg model integral: " << bkgModel->createIntegral(RooArgList(*mass,*cosT))->getVal() <<std::endl;
    //build combined fit model
    RooAddPdf *comb = new RooAddPdf(Form("Data_%s_TEMPLATE2DFIT_%s",mcName.Data(),catIt->Data()),"",RooArgList(*signalHistPdf,*bkgModel),
				    RooArgList(*thisNsig,*thisNbkg) );

    
    pdfsToDelete.push_back(bkgModel);
    pdfsToDelete.push_back(comb);

    combFit->addPdf(*comb,*catIt);
    delete signalData;
    //delete signalDataHist;    
  }    

  RooDataSet *ds = (RooDataSet*)ws->data("Data_Combined");

  combFit->fitTo(*ds,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=combFit->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_%s_TEMPLATE2DFIT_fitResult",mcName.Data()) );
  
  ws->import(*combFit);
  ws->import(*res);

  while(!pdfsToDelete.empty()){ delete pdfsToDelete.back();  pdfsToDelete.pop_back(); }
  delete combFit;
}


void MakeSpinFits::Make2DFloatingSignalTest(TString massMcName,TString costMcName){
  std::cout << "MakeSpinFits::Make2DFloatingSignalTest" <<std::endl;
  if(ws==0) return;

  RooCategory *cat = ((RooCategory*)ws->obj("evtcat"));

  RooSimultaneous *combFit = new RooSimultaneous(Form("Data_m_%s_c_%s_IND2DFIT",massMcName.Data(),costMcName.Data()),"",*cat);

  RooRealVar *nSig = new RooRealVar( Form("Data_m_%s_c_%s_IND2DFIT_Nsig",massMcName.Data(),costMcName.Data() ), "", 0, 0,1e9);
  RooRealVar *nBkg = new RooRealVar( Form("Data_m_%s_c_%s_IND2DFIT_Nbkg",massMcName.Data(),costMcName.Data() ), "", 0, 0,1e9);

  RooRealVar *mass = ws->var("mass");
  RooRealVar *cosT = ws->var("cosT");



  std::vector<TString>::const_iterator catIt = catLabels.begin();
  for(; catIt != catLabels.end(); catIt++){
    cout << *catIt <<endl;

    RooAbsPdf *signalMassModel = ws->pdf( Form("Data_%s_FIT_%s",massMcName.Data(),catIt->Data()));
    RooAbsPdf *signalCosModel = ws->pdf(Form("%s_FIT_%s_cosTpdf",costMcName.Data(),catIt->Data()));
    
    RooAbsPdf *bkgMassModel    = ws->pdf( Form("Data_BKGFIT_%s_bkgShape",catIt->Data() ) );

    RooAbsPdf *bkgCosModel    = ws->pdf(Form("Data_BKGFIT_%s_cosTpdf",catIt->Data()));

    double totalEventsB = ws->data("Data_Combined")->sumEntries();
    double thisCatB     = ws->data("Data_Combined")->sumEntries(TString("evtcat==evtcat::")+*catIt);

    RooRealVar *nBkg = new RooRealVar( Form("Data_m_%s_c_%s_IND2DFIT_%s_Nbkg",massMcName.Data(),costMcName.Data(),catIt->Data()), "", thisCatB,0,1e9 ); 
    //float signal yield
    RooRealVar *nSig    = new RooRealVar( Form("Data_m_%s_c_%s_IND2DFIT_%s_Nsig",massMcName.Data(),costMcName.Data(),catIt->Data() ), "", 0,0,1e9);

    //RooExtendPdf *exSignalModel = new RooExtendPdf(Form("Data_%s_IND2DFIT_%s_signalModel",mcName.Data(),catIt->Data()),"",*signalModel,*thisNsig);
    /*
    std::cout << "sig model MASS integral: " << signalMassModel->createIntegral(*mass)->getVal() <<std::endl;
    std::cout << "bkg model MASS integral: " << bkgMassModel->createIntegral(*mass)->getVal() <<std::endl;

    std::cout << "sig model COST integral: " << signalCosModel->createIntegral(*cosT)->getVal() <<std::endl;
    std::cout << "bkg model COST integral: " << bkgCosModel->createIntegral(*cosT)->getVal() <<std::endl;
    */
    RooProdPdf *bkgModel = new RooProdPdf(Form("Data_m_%s_c_%s_IND2DFIT_%s_bkgModel",massMcName.Data(),costMcName.Data(),catIt->Data()),"",
					     RooArgList(*bkgMassModel,*bkgCosModel));
    
    RooProdPdf *signalModel = new RooProdPdf(Form("Data_m_%s_c_%s_IND2DFIT_%s_signalModel",massMcName.Data(),costMcName.Data(),catIt->Data()),"",
					     RooArgList(*signalMassModel,*signalCosModel));
    /*
    std::cout << "sig model integral: " << signalModel->createIntegral(RooArgList(*mass,*cosT))->getVal() <<std::endl;
    std::cout << "bkg model integral: " << bkgModel->createIntegral(RooArgList(*mass,*cosT))->getVal() <<std::endl;
    */
    RooAddPdf *comb = new RooAddPdf(Form("Data_m_%s_c_%s_IND2DFIT_%s",massMcName.Data(),costMcName.Data(),catIt->Data()),"",RooArgList(*signalModel,*bkgModel),
				    RooArgList(*nSig,*nBkg) );

    combFit->addPdf(*comb,*catIt);
  }    

  RooDataSet *ds = (RooDataSet*)ws->data("Data_Combined");

  combFit->fitTo(*ds,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=combFit->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_m_%s_c_%s_IND2DFIT_fitResult",massMcName.Data(),costMcName.Data()) );
  
  ws->import(*combFit);
  ws->import(*res);

}
/*
RooAbsPdf* MakeSpinFits::Make2DSignalModel(TString massMcName,TString costMcName,TString catTag,TString inset){
  return 0;
  /*
  RooAbsPdf *signalMassModel = ws->pdf( Form("Data_%s_FIT_%s",massMcName.Data(),catTag.Data()));
  RooAbsPdf *signalCosModel = ws->pdf(Form("%s_FIT_%s_cosTpdf",costMcName.Data(),catTag.Data()));
    TString name = TString("Data_m_")+massMcName+"_c_"+costMcName+"_"+inset+"_"+catTag+"_signalModel";
    if(!use2DIntegralMorph) return new RooProdPdf(name,"",RooArgList(*signalMassModel,*signalCosModel));
    
    RooCategory *cosThetaBins = (RooCategory*)ws->obj("CosThetaBins");
    
    cosThetaBins->setIndex(0);
    TString low = cosThetaBins->getLabel();
    cosThetaBins->setIndex(cosThetaBins->numBins("")-1);
    TString high = cosThetaBins->getLabel();
    

    RooAbsPdf *signalMassModel_low = ws->pdf( TString("Data_")+massMcName+"_FIT_"+catTag+"_"+low );
    RooAbsPdf *signalMassModel_high = ws->pdf( TString("Data_")+massMcName+"_FIT_"+catTag+"_"+high );

    RooRealVar *mass = ws->var("mass");
    RooRealVar *cosT = ws->var("cosT");

    mass->setBins(400,"cache");
    cosT->setBins(100,"cache");

    RooIntegralMorph *morph = new RooIntegralMorph(name,"",*mass,*cosT);
    morph->setCacheAlpha(kTRUE);
    return morph;

}

RooAbsPdf* MakeSpinFits::Make2DBkgModel(TString massMcName,TString costMcName,TString catTag,TString inset){return 0;}
*/
void MakeSpinFits::MakeBackgroundOnlyFit(TString catTag, float cosTlow, float cosThigh){
  std::cout << "MakeSpinFits::MakeBackgroundOnlyFit" <<std::endl;
  if(ws==0) return;
  RooRealVar mass = *(ws->var("mass"));
  cout << "1" <<endl;

  RooAbsData *ds = ws->data("Data_Combined")->reduce(TString("evtcat==evtcat::")+catTag);

  TString outputTag = catTag;
  if(cosTlow > -1 || cosThigh < 1){
    outputTag+=Form("_cosT_%0.2f_%0.2f",cosTlow,cosThigh);
    ds = (RooDataSet*)ds->reduce(  Form("cosT >= %0.3f && cosT < %0.3f",cosTlow,cosThigh) );
  }


  //background model
  RooAbsPdf* BkgShape;
 
  switch(fitType){
  case kExp:
    {
      //double exponential
    RooRealVar* alpha1 = new RooRealVar(Form("Data_BKGFIT_%s_alpha1",outputTag.Data()),"alpha1",-0.1,-1.,0.);
    RooRealVar* alpha2 = new RooRealVar(Form("Data_BKGFIT_%s_alpha2",outputTag.Data()),"alpha2",-0.1,-1.,0.);
    RooRealVar* f_bkg  = new RooRealVar(Form("Data_BKGFIT_%s_f",outputTag.Data()),"f_bkg",0.1,0,1);
    RooExponential* exp1 = new RooExponential(Form("Data_BKGFIT_%s_exp1",outputTag.Data()),"exp1",mass,*alpha1);
    RooExponential* exp2 = new RooExponential(Form("Data_BKGFIT_%s_exp2",outputTag.Data()),"exp2",mass,*alpha2);
    
    BkgShape = new RooAddPdf(Form("Data_BKGFIT_%s_bkgShape",outputTag.Data()),"Background Model",
			     RooArgList(*exp1,*exp2),*f_bkg);
    break;
    }
  case kPoly:
    {
      //5th order polynomial
    RooRealVar *pC = new RooRealVar(Form("Data_BKGFIT_%s_pC",outputTag.Data()),"pC",1);
    RooRealVar *p0 = new RooRealVar(Form("Data_BKGFIT_%s_p0",outputTag.Data()),"p0",0,-10,10);
    RooRealVar *p1 = new RooRealVar(Form("Data_BKGFIT_%s_p1",outputTag.Data()),"p1",0,-10,10);
    RooRealVar *p2 = new RooRealVar(Form("Data_BKGFIT_%s_p2",outputTag.Data()),"p2",0,-10,10);
    RooRealVar *p3 = new RooRealVar(Form("Data_BKGFIT_%s_p3",outputTag.Data()),"p3",0,-10,10);
    RooRealVar *p4 = new RooRealVar(Form("Data_BKGFIT_%s_p4",outputTag.Data()),"p4",0,-10,10);
    //enforce all coefficients positive
    RooFormulaVar *pCmod = new RooFormulaVar(Form("Data_BKGFIT_%s_pCmod",outputTag.Data()),"","@0*@0",*pC);
    RooFormulaVar *p0mod = new RooFormulaVar(Form("Data_BKGFIT_%s_p0mod",outputTag.Data()),"","@0*@0",*p0);
    RooFormulaVar *p1mod = new RooFormulaVar(Form("Data_BKGFIT_%s_p1mod",outputTag.Data()),"","@0*@0",*p1);
    RooFormulaVar *p2mod = new RooFormulaVar(Form("Data_BKGFIT_%s_p2mod",outputTag.Data()),"","@0*@0",*p2);
    RooFormulaVar *p3mod = new RooFormulaVar(Form("Data_BKGFIT_%s_p3mod",outputTag.Data()),"","@0*@0",*p3);
    RooFormulaVar *p4mod = new RooFormulaVar(Form("Data_BKGFIT_%s_p4mod",outputTag.Data()),"","@0*@0",*p4);

  

    RooArgList *args;
    if(cosTlow > -1 || cosThigh < 1){
      args = new RooArgList(*pCmod,*p0mod,*p1mod,*p2mod,*p3mod);
    }else{
      args = new RooArgList(*pCmod,*p0mod,*p1mod,*p2mod,*p3mod,*p4mod);
    }

    BkgShape = new RooBernstein(Form("Data_BKGFIT_%s_bkgShape",outputTag.Data()),"Background Model",mass,*args);
    break;
    }
  default:
    std::cout << "INVALID BACKGROUND MODEL" << std::endl;
    assert(false);
    break;
  }

  RooRealVar *Nbkg = new RooRealVar(Form("Data_BKGFIT_%s_Nbkg",outputTag.Data()),"N background Events",ds->sumEntries(),0,1e9);
  //extended fit model
  RooExtendPdf *BkgModel = new RooExtendPdf(Form("Data_BKGFIT_%s_bkgModel",outputTag.Data()),"Background Model",*BkgShape,*Nbkg);
  
  //  mass.setRange("fitLow", 100,120.5);
  //mass.setRange("fitHigh",127.5,180);

  BkgModel->fitTo(*ds,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=BkgModel->fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_BKGFIT_%s_fitResult",outputTag.Data()) );

  //set all parameters of the background fit as constant after the fit
  RooArgSet* vars = BkgModel->getVariables();
  RooFIter iter = vars->fwdIterator();
  RooAbsArg* a;
  while( (a = iter.next()) ){
    if(a->GetName() == "mass") continue;
    static_cast<RooRealVar*>(a)->setConstant(kTRUE);
  }

  //make cosT pdf
  if(!(cosTlow > -1 || cosThigh < 1)){
    RooRealVar* cosT = ws->var("cosT");
    RooDataSet *dsred = (RooDataSet*)ds->reduce("mass<120 || mass > 130");
    RooDataHist hist(Form("Data_BKGFIT_%s_cosThist",outputTag.Data()),"Data Hist for cos(theta)",RooArgSet(*cosT),*ds);
    RooHistPdf cosTkde(Form("Data_BKGFIT_%s_cosTpdf",outputTag.Data()),"Hist PDF for cos(theta)",RooArgSet(*cosT),hist);
    
    ws->import(cosTkde);
  }

  ws->import(*BkgModel);
  ws->import(*Nbkg);
  ws->import(*res);

  delete ds;
}

void MakeSpinFits::getSimpleBkgSubtraction(TString mcName,TString tag){
  std::cout << "\n\n\nMakeSpinFits::getSimpleBkgSubtraction\n\n\n" <<std::endl;
  RooDataSet *data = (RooDataSet*)ws->data("Data_Combined")->reduce(TString("evtcat==evtcat::")+tag);
  RooRealVar *mass = ws->var("mass");
  RooRealVar *cosT = ws->var("cosT");

  cosT->setBins(10);

  double se  = ws->var( Form("%s_FIT_%s_sigmaEff",mcName.Data(),tag.Data()))->getVal();
  double m   = ws->var( Form("Data_%s_FIT_mean",mcName.Data()))->getVal();

  double nSig = ((RooFormulaVar*)ws->obj( Form("Data_%s_FULLFIT_%s_Nsig",mcName.Data(),tag.Data()) ))->getVal()*0.682; //only take 68% of the signal, since we only do +- 1 sigma_eff


  /*
  RooAbsPdf *bkgPdf = ws->pdf(Form("Data_BKGFIT_%s_bkgModel",tag.Data()));
  double nBkgTot = ((RooRealVar*)ws->obj(Form("Data_%s_FULLFIT_%s_Nbkg",mcName.Data(),tag.Data())))->getVal();
  RooRealVar range("range","",m-se,m+se);
  RooRealVar all("all","",100,170);
  double nBkg  = bkgPdf->createIntegral(range)->getVal()/bkgPdf->createIntegral(all)->getVal()*nBkgTot;
  */

  RooDataSet *bkg = (RooDataSet*)data->reduce( Form("mass < %f || mass > %f",m-3*se,m+3*se) );
  RooDataHist bkgCos("bkgCos","",RooArgSet(*cosT),*bkg);
  RooHistPdf  bkgCosPdf("bkgCosPdf","",RooArgSet(*cosT),bkgCos);

  RooDataSet *sig = (RooDataSet*)data->reduce( Form("mass < %f && mass > %f",m+se,m-se) );
  RooDataHist sigHist("sigCos","",RooArgSet(*cosT),*sig);

  double nBkg = sig->sumEntries()-nSig;

  std::cout << m << "    " << se << "   " << nBkg << "   " << sig->sumEntries()<< std::endl;
  
  RooDataHist *sigBkgSub = subtract(*cosT,sigHist,bkgCosPdf,nBkg);
  sigBkgSub->SetName(Form("Data_%s_%s_bkgSub_cosT",mcName.Data(),tag.Data()));

  ws->import(*sigBkgSub);
  delete data;
}

void MakeSpinFits::getSimpleTotalBkgSubtraction(TString mcName){
  std::cout << "\n\n\nMakeSpinFits::getSimpleTotalBkgSubtraction\n\n\n" <<std::endl;
  RooDataSet *data = (RooDataSet*)ws->data( "Data_Combined" );
  RooRealVar *mass = ws->var("mass");
  RooRealVar *cosT = ws->var("cosT");

  cosT->setBins(10);

  double se  = ws->var( Form("%s_FIT_Combined_sigmaEff",mcName.Data()))->getVal();
  double m   =  ws->var( Form("Data_%s_FIT_mean",mcName.Data()))->getVal();

  double nSig = ws->var( Form("Data_%s_FULLFIT_Nsig",mcName.Data()) )->getVal();
  
  
  RooDataSet *bkg = (RooDataSet*)data->reduce( Form("mass < %f || mass > %f",m-3*se,m+3*se) );
  RooDataHist bkgCos("bkgCos","",RooArgSet(*cosT),*bkg);
  RooHistPdf  bkgCosPdf("bkgCosPdf","",RooArgSet(*cosT),bkgCos);

  RooDataSet *sig = (RooDataSet*)data->reduce( Form("mass < %f && mass > %f",m+se,m-se) );
  RooDataHist sigHist("sigCos","",RooArgSet(*cosT),*sig);

  double nBkg = sig->sumEntries()-nSig;

  std::cout << m << "    " << se << "   " << nBkg << "   " << sig->sumEntries()<< std::endl;
  
  RooDataHist *sigBkgSub = subtract(*cosT,sigHist,bkgCosPdf,nBkg);
  sigBkgSub->SetName(Form("Data_%s_Combined_bkgSub_cosT",mcName.Data()));

  ws->import(*sigBkgSub);
}

void MakeSpinFits::run(){
  std::cout << "\n\n\nMakeSpinFits::run\n\n\n" <<std::endl;
  if(ws==0) return;
  RooMsgService::instance().setSilentMode(true);

  RooCategory cosThetaBins("CosThetaBins","CosThetaBins");
  for(int i=0;i<NcosTbins;i++){
    cosThetaBins.defineType(Form("cosT_%0.2f_%0.2f",cosTbinEdges[i],cosTbinEdges[i+1]),cosThetaBins.numBins(""));
  }
  ws->import(cosThetaBins);

  binDatasetCosT(*(ws->data("Data_Combined")),"Data");

  //run fits in the correct order for each MC type
  std::vector<TString>::const_iterator mcIt = mcLabel.begin();
  for(; mcIt != mcLabel.end(); mcIt++){
    binDatasetCosT(*(ws->data(*mcIt+"_Combined")),*mcIt);

    std::vector<TString>::const_iterator catIt = catLabels.begin();
    MakeSignalFit("Combined",*mcIt);
    for(int i=0;i<NcosTbins;i++) MakeSignalFit("Combined",*mcIt,cosTbinEdges[i],cosTbinEdges[i+1]);

    MakeSignalFitForFit("Combined",*mcIt);

    for(; catIt != catLabels.end(); catIt++){
      MakeSignalFit(*catIt,*mcIt);
      MakeSignalFitForFit(*catIt,*mcIt);
      for(int i=0;i<NcosTbins;i++){
	MakeSignalFit(*catIt,*mcIt,cosTbinEdges[i],cosTbinEdges[i+1]);
	MakeSignalFitForFit(*catIt+Form("_cosT_%0.2f_%0.2f",cosTbinEdges[i],cosTbinEdges[i+1]),*mcIt);
      }
      if(mcIt == mcLabel.begin()){
	MakeBackgroundOnlyFit(*catIt);
	for(int i=0;i<NcosTbins;i++) MakeBackgroundOnlyFit(*catIt,cosTbinEdges[i],cosTbinEdges[i+1]);	
      }
    }

    MakeCombinedSignalTest(*mcIt);
    MakeFloatingSignalTest(*mcIt);

    MakeCombinedSignalTest(*mcIt,true);

    catIt = catLabels.begin();
    for(; catIt != catLabels.end(); catIt++){
      getSimpleBkgSubtraction(*mcIt,*catIt);
    }
    MakeCombinedSignalSpin(*mcIt);
    getSimpleTotalBkgSubtraction(*mcIt);
    AddCombinedBkgOnlySWeight(*mcIt);
    std::cout << "DONE WITH " << *mcIt <<std::endl;
    ws->Write(ws->GetName(),TObject::kWriteDelete);
  }
    std::cout << "DONE" <<std::endl;

  //make the 2D fits for all possible combinations of cos(theta) pdf and lineshape
  /*
  for(int i=0;i<mcLabel.size();i++){
    for(int j=0;j<mcLabel.size();j++){
      Make2DCombinedSignalTest(mcLabel.at(i),mcLabel.at(j));
      Make2DFloatingSignalTest(mcLabel.at(i),mcLabel.at(j));
    }
  }
  */ // combinatorics too high for more than 2-3 signal samples
  
}

int MakeSpinFits::specifySamples(std::vector<std::string> samples){
  mcLabel.clear();

  RooCategory* labels = ((RooCategory*)ws->obj("labels"));
  
  for(std::vector<std::string>::const_iterator it = samples.begin();
      it != samples.end(); it++)
    {
      std::cout << "Running: " << *it << std::endl;
      if( labels->setLabel(it->c_str(),kFALSE) != 0) //MC sample not in workspace
	{
	  std::cout << "\n\nERROR: MC Sample " << *it 
		    << " not in workspace" <<std::endl;
	  return -1;
	}
      mcLabel.push_back(TString(it->c_str()));
    }
  return 0;
}

void MakeSpinFits::save(){
  if(ws==0 || outputFile==0) return;
  std::cout << "SAVING" <<std::endl;
  outputFile->cd();
  std::cout <<"WRITING" <<std::endl;
  RooCategory *outLabels = new RooCategory("fitlabels","");
  for(std::vector<TString>::const_iterator it = mcLabel.begin();
      it != mcLabel.end(); it++){
    outLabels->defineType(it->Data(),it-mcLabel.begin());
  }
  ws->import(*outLabels);
  ws->Write(ws->GetName(),TObject::kWriteDelete);
  std::cout << "CLOSING" <<std::endl;
  outputFile->Close();
}

void MakeSpinFits::binDatasetCosT(RooAbsData& data,TString name){
  RooCategory *cosTBins = (RooCategory*)ws->obj("CosThetaBins");
  RooCategory *baseCats = (RooCategory*)ws->obj("evtcat");

  RooRealVar *mass = ws->var("mass");
  RooRealVar *cosT = ws->var("cosT");

  cosTcatLabels.clear();

  RooCategory *cat = new RooCategory("evtcat_cosT","evtcat_cosT");
  for(int iCat=0;iCat<baseCats->numBins("");iCat++){
    baseCats->setIndex(iCat);
    TString catName = baseCats->getLabel();
    for(int iCos=0;iCos<cosTBins->numBins("");iCos++){
      cosTBins->setIndex(iCos);
      TString cosName = cosTBins->getLabel();
      cosTcatLabels.push_back(catName+"_"+cosName);
      cat->defineType(catName+"_"+cosName, cat->numBins(""));
    }
  }
  RooDataSet* dataComb = new RooDataSet(name+"_Combined_CosTBin","",RooArgSet(*mass,*cosT,*cat) );
  for(int iCat=0;iCat<baseCats->numBins("");iCat++){
    baseCats->setIndex(iCat);
    TString catName = baseCats->getLabel();
    for(int iCos=0;iCos<cosTBins->numBins("");iCos++){
      cosTBins->setIndex(iCos);
      TString cosName = cosTBins->getLabel();
      std::pair<float,float> cosRange = getCosTRangeFromCatName(cosName);

      RooDataSet * d = (RooDataSet*)data.reduce( Form("evtcat==%d && cosT >= %0.3f && cosT < %0.3f",iCat,cosRange.first,cosRange.second) );
      RooDataSet *tmp = new RooDataSet("DataCat"+catName+"_"+cosName,"",RooArgSet(*mass,*cosT),RooFit::Index(*cat),RooFit::Import(Form("%s_%s",catName.Data(),cosName.Data()),*d));
      dataComb->append(*tmp);
      delete tmp;
      delete d;
    }
  }
  ws->import(*cat);
  ws->import(*dataComb);
  delete dataComb;
  delete cat;
}

void MakeSpinFits::setCosTBins(const int N, const float *edges){
  if(N > 100) assert(false);
  if(edges==0) assert(false);
  //if(cosTbinEdges) delete cosTbinEdges;
  cosTbinEdges = new float[N];
  for(int i=0;i<N;i++){
    cosTbinEdges[i] = edges[i];
  }
  NcosTbins=N-1;
}

std::pair<float,float> MakeSpinFits::getCosTRangeFromCatName(TString name){
  TObjArray *tokens = name.Tokenize("_");
  float low = atof(tokens->At(1)->GetName());
  float high = atof(tokens->At(2)->GetName());
  delete tokens;
  return std::pair<float,float> (low,high);
}
