#include "MakeSpinToy.h"
#include <iostream>

MakeSpinToy::MakeSpinToy(TString fileName,TString wsName)
  //nCat(2);
{
  TFile *f = new TFile(fileName);
  ws = (RooWorkspace*)f->Get(wsName);

  MakeSpinFits::getLabels("evtcat",&catLabels,ws);
  MakeSpinFits::getLabels("evtcat_cosT",&catCosTLabels,ws);
  
  fitType=MakeSpinFits::kPoly;

  mass = ws->var("mass");
  cosT = ws->var("cosT");
  //cosT->setRange(0,1);
  cosT->setBins(10);
  S = new RooRealVar("S","",0,-1e7,1e7);

  GenMinusFit = new RooRealVar("NgenMinusNfit","",0,-1e6,1e6);

  S_TruthHgg = new RooDataSet("S_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_TruthALT = new RooDataSet("S_TruthALT","",RooArgSet(*S,*GenMinusFit));

  S_splot_TruthHgg = new RooDataSet("S_splot_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_splot_TruthALT = new RooDataSet("S_splot_TruthALT","",RooArgSet(*S,*GenMinusFit));

  S_2D_TruthHgg = new RooDataSet("S_2D_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_2D_TruthALT = new RooDataSet("S_2D_TruthALT","",RooArgSet(*S,*GenMinusFit));

  S_tot_TruthHgg = new RooDataSet("S_tot_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_tot_TruthALT = new RooDataSet("S_tot_TruthALT","",RooArgSet(*S,*GenMinusFit));

  S_2DFIT_TruthHgg = new RooDataSet("S_2DFIT_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_2DFIT_TruthALT = new RooDataSet("S_2DFIT_TruthALT","",RooArgSet(*S,*GenMinusFit));

  S_2DTEMPFIT_TruthHgg = new RooDataSet("S_2DTEMPFIT_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_2DTEMPFIT_TruthALT = new RooDataSet("S_2DTEMPFIT_TruthALT","",RooArgSet(*S,*GenMinusFit));

  S_2DBINFIT_TruthHgg = new RooDataSet("S_2DBINFIT_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_2DBINFIT_TruthALT = new RooDataSet("S_2DBINFIT_TruthALT","",RooArgSet(*S,*GenMinusFit));

  S_TruthData=0;
  S_tot_TruthData=0;
  S_splot_TruthData=0;
  S_2D_TruthData=0;
  S_2DFIT_TruthData=0;
  S_2DTEMPFIT_TruthData=0;
  S_2DBINFIT_TruthData=0;

  setSaveWorkspaces(false);
  setDoData(false);

  mcLabels[0]="Data";
  mcLabels[1]="Hgg125";
  mcLabels[2]="RSG125";
}

void MakeSpinToy::setDoData(bool b){
  doData=b;
  if(doData){
    S_TruthData = new RooDataSet("S_TruthData","",RooArgSet(*S,*GenMinusFit));
    S_splot_TruthData = new RooDataSet("S_splot_TruthData","",RooArgSet(*S,*GenMinusFit));
    S_2D_TruthData = new RooDataSet("S_2D_TruthData","",RooArgSet(*S,*GenMinusFit));
    S_2DFIT_TruthData = new RooDataSet("S_2DFIT_TruthData","",RooArgSet(*S,*GenMinusFit));
    S_2DTEMPFIT_TruthData = new RooDataSet("S_2DFIT_TruthData","",RooArgSet(*S,*GenMinusFit));
    S_2DBINFIT_TruthData = new RooDataSet("S_2DFIT_TruthData","",RooArgSet(*S,*GenMinusFit));
    S_tot_TruthData = new RooDataSet("S_tot_TruthData","",RooArgSet(*S,*GenMinusFit));
    
  }
}

void MakeSpinToy::generateToyWorkspace(RooWorkspace *toyws,genType gen){
  RooCategory *cat = new RooCategory("evtcat","evtcat");
  int i=0;
  for( std::vector<TString>::const_iterator it = catLabels.begin();
       it != catLabels.end(); it++,i++){
    std::cout << *it << "   " << i << std::endl;
    cat->defineType( *it,i);
  }
  RooDataSet* dataComb = new RooDataSet("Data_Combined","",RooArgSet(*mass,*cosT,*cat) );

  TRandom3 rng(0);
  float Nsig = rng.Poisson( ws->data(mcLabels[1]+"_Combined")->sumEntries()*targetLumi/nominalLumi); 
  //float Nsig = rng.Poisson( 607.*targetLumi/12.);
  for( std::vector<TString>::const_iterator it = catLabels.begin();
       it != catLabels.end(); it++,i++){
    std::cout << *it <<std::endl;
    generateToyWorkspace(toyws,it->Data(),gen,Nsig);
    RooDataSet* d = (RooDataSet*)toyws->data( Form("Data_%s",it->Data()) );
    RooDataSet* tmp = new RooDataSet("DataCat+"+*it,"",RooArgSet(*mass,*cosT),RooFit::Index(*cat),RooFit::Import(*it,*d) );
    dataComb->append(*tmp);
    std::cout << Form("%s_FIT_%s_sigmaEff",mcLabels[1].Data(),it->Data()) <<"     " <<ws << std::endl;

    toyws->import(*(ws->var(Form("%s_FIT_%s_sigmaEff",mcLabels[1].Data(),it->Data()) )));
    toyws->import(*(ws->var(Form("%s_FIT_%s_sigmaEff",mcLabels[2].Data(),it->Data()) )));
    delete tmp;
  }
  toyws->import(*cat);
  toyws->import(*dataComb);
  toyws->import(*(ws->data(mcLabels[1]+"_Combined")));
  toyws->import(*(ws->data(mcLabels[2]+"_Combined")));
  toyws->import(*(ws->var(mcLabels[1]+"_FIT_Combined_sigmaEff")));
  toyws->import(*(ws->var(mcLabels[2]+"_FIT_Combined_sigmaEff")));
  toyws->import(*((RooCategory*)ws->obj("labels")));

  delete dataComb;
}

void MakeSpinToy::generateToyWorkspace(RooWorkspace *toyws,const char* cat,genType gen,float nSigTot){
  std::cout << "generateToyWorkspace per cat" <<std::endl;
  //get the generation PDFs
  RooCategory* cosTCats = (RooCategory*)ws->obj("CosThetaBins");

  RooDataSet *data = new RooDataSet(Form("Data_%s",cat),"",RooArgSet(*mass,*cosT));
  RooDataSet *toysig = new RooDataSet(Form("ToySigData_%s",cat),"",RooArgSet(*mass,*cosT));
  int sigGen=0,bkgGen=0;
  RooAbsPdf* hggMassPdf = ws->pdf( Form("%s_FIT_%s",mcLabels[1].Data(),cat) );
  RooAbsPdf* altMassPdf = ws->pdf( Form("%s_FIT_%s",mcLabels[2].Data(),cat) );


  for(int iCos=0;iCos<cosTCats->numBins("");iCos++){
    cosTCats->setIndex(iCos);
    const char* csCat = Form("%s_%s",cat,cosTCats->getLabel());
    //mass

    //ws->var(Form("%s_FIT_%s_mean",mcLabels[1].Data(),cat))->setVal(125.);
    //ws->var(Form("%s_FIT_%s_mean",mcLabels[2].Data(),cat))->setVal(125.);

    RooAbsPdf* bkgMassPdf = ws->pdf( Form("Data_BKGFIT_%s_bkgModel",csCat) );
  
    cosT->setBins(50,"generate");
    mass->setBins(200,"generate");
    //cosT
    RooDataSet *tmp  = (RooDataSet*)(ws->data("Data_Combined_CosTBin")->reduce(Form("evtcat_cosT==evtcat_cosT::%s && ((mass>110 && mass<120) || (mass>130 && mass<140))",csCat)));
    RooDataSet *tmpH = (RooDataSet*)ws->data(Form("%s_Combined_CosTBin",mcLabels[1].Data()))->reduce(Form("evtcat_cosT==evtcat_cosT::%s",csCat));
    RooDataSet *tmpR = (RooDataSet*)ws->data(Form("%s_Combined_CosTBin",mcLabels[2].Data()))->reduce(Form("evtcat_cosT==evtcat_cosT::%s",csCat));

    RooDataHist histBkg ("bkgDataHist","",*cosT,"generate"); histBkg.add(*tmp);
    RooDataHist histHgg ("hggDataHist","",RooArgSet(*mass,*cosT),"generate"); histHgg.add(*tmpH);
    RooDataHist histAlt ("altDataHist","",RooArgSet(*mass,*cosT),"generate"); histAlt.add(*tmpR);
    
    RooHistPdf bkgGenPdf("bkgGenPdf","",*cosT,histBkg);
    RooHistPdf hggGenPdf("hggGenPdf","",RooArgSet(*mass,*cosT),histHgg);
    RooHistPdf altGenPdf("altGenPdf","",RooArgSet(*mass,*cosT),histAlt);
    
    /*
      RooDataHist histBkg ("bkgDataHist","",*cosT,"generate"); histBkg.add(*tmp);
      RooDataHist histHgg ("hggDataHist","",*cosT,"generate"); histHgg.add(*tmpH);
      RooDataHist histAlt ("altDataHist","",*cosT,"generate"); histAlt.add(*tmpR);
      
      RooHistPdf bkgGenPdf("bkgGenPdf","",*cosT,histBkg);
      RooHistPdf hggGenPdf("hggGenPdf","",*cosT,histHgg);
      RooHistPdf altGenPdf("altGenPdf","",*cosT,histAlt);
    */
    /*
      RooKeysPdf* bkgGenPdf = new RooKeysPdf("bkgGenPdf","",*cosT,*tmp);
      RooKeysPdf* hggGenPdf = new RooKeysPdf("hggGenPdf","",*cosT,*tmpH);
      RooKeysPdf* altGenPdf = new RooKeysPdf("altGenPdf","",*cosT,*tmpR);
    */
    /*
      tmp->SetName(Form("OriginalData_%s",cat));
      RooDataHist tmpHist("bkgGenHist","",*cosT,*tmp);
      RooHistPdf* bkgGenPdf = new RooHistPdf("bkgGenPdf","",*cosT,tmpHist);
      
      std::cout << "BKG GEN PDF:  " << bkgGenPdf <<std::endl;
      if(!bkgGenPdf) return;
      RooHistPdf* hggGenPdf = (RooHistPdf*)ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",cat));
      RooHistPdf* altGenPdf = (RooHistPdf*)ws->pdf(Form("ALT125_FIT_%s_cosTpdf",cat));
    */


    //get the expected number of events in this categoy
    //Data_Hgg125_FULAALFIT_EB_0_Nbkg
    int Nbkg = ((RooFormulaVar*)ws->obj(Form("Data_%s_FULLCOSTFIT_%s_Nbkg",mcLabels[gen].Data(),csCat)))->getVal()*targetLumi/nominalLumi;
    int Nsig = nSigTot * ws->var(Form("Data_%s_FULLCOSTFIT_%s_fsig",mcLabels[gen].Data(),csCat))->getVal();

    TRandom3 rng(0);
    int thisNbkg = (int)rng.Poisson(Nbkg);
    int thisNsig = Nsig;//(int)rng.Poisson(Nsig); //number of events to generate
    
    std::cout << "\nGenerating  " << thisNbkg << " background and " << thisNsig << " signal events\n" <<std::endl;

    sigGen+=thisNsig;
    bkgGen+=thisNbkg;


    RooDataSet *M = bkgMassPdf->generate(*mass,thisNbkg);
    RooDataSet *C = bkgGenPdf.generate(*cosT,thisNbkg,RooFit::AutoBinned(kFALSE));
    
    
    RooDataSet *mcM, *mcC;
    RooDataSet *mc;
    switch(gen){
    case Hgg125:
      std::cout << "GENERATING HIGGS SPIN DISTRIBUTION" <<std::endl;
      //mcM = hggMassPdf->generate(*mass,thisNsig);
      //mcC = hggGenPdf.generate(*cosT,thisNsig,RooFit::AutoBinned(kFALSE));
      mc = hggGenPdf.generate(RooArgSet(*mass,*cosT),thisNsig,RooFit::AutoBinned(kFALSE));
      break;
      
    case ALT125:
      std::cout << "GENERATING ALTERNATE SPIN DISTRIBUTION" <<std::endl;
      //mcM = altMassPdf->generate(*mass,thisNsig); 
      //mcC = altGenPdf.generate(*cosT,thisNsig,RooFit::AutoBinned(kFALSE));
      mc = altGenPdf.generate(RooArgSet(*mass,*cosT),thisNsig,RooFit::AutoBinned(kFALSE));
      break;
      
    default:
      return;
    }
    std::cout << "M size: " << M->sumEntries() << std::endl << "C size: " << C->sumEntries() <<std::endl;
    //std::cout << "M size: " << mcM->sumEntries() << std::endl << "C size: " << mcC->sumEntries() <<std::endl;
    
    M->Print();
    C->Print();
    //mcM->Print();
    //mcC->Print();
    mc->Print();
    M->merge(C);
    //mcM->merge(mcC);
    M->append(*mc);
    
    data->append(*M);
    toysig->append(*mc);


    delete tmp;
    delete tmpH;
    delete tmpR;
  }  

  RooRealVar nSigGen(Form("N_gen_sig_%s",cat),"",sigGen);
  RooRealVar nBkgGen(Form("N_gen_bkg_%s",cat),"",bkgGen);

  toyws->import(*data);
  toyws->import(*toysig);
  toyws->import(*hggMassPdf);
  toyws->import(*altMassPdf);
  toyws->import(nSigGen);
  toyws->import(nBkgGen);

  delete data;
  delete toysig;

  //delete hggMassPdf;
  //delete altMassPdf;

  std::cout << "Done Generating Category " << cat << std::endl;
}
 
double MakeSpinToy::computeLL(RooAbsPdf* pdf, RooAbsData* data,RooRealVar* var, int rebin){
  //TH1F dataHist("dataHist","",nBins,-1,1);

  //data->fillHistogram(&dataHist,*var);

  TH1F* dataHist = getHistogram(data,"dataHist",rebin);
  TH1F* pdfHist = (TH1F*)pdf->createHistogram("cosT",var->getBins()/rebin);

  std::cout << dataHist->Integral() << std::endl;
  pdfHist->Scale(dataHist->Integral()/pdfHist->Integral());


  std::cout << "Bins:" <<std::endl;
  for(int i=0;i<cosT->getBins();i++){
    std::cout << "\t" << i+1 << ":    " << dataHist->GetBinContent(i+1) << " +- " << dataHist->GetBinError(i+1) << "          pdf: " 
	      << pdfHist->GetBinContent(i+1) << std::endl;
  }

  double p = dataHist->Chi2Test(pdfHist,"WW");

  std::cout << "p: " << p << std::endl;

  delete pdfHist;
  delete dataHist;
  if(p==0) return 0;
  return TMath::Log(p);

  /*
  RooDataHist dataHist("tmp","",*var,*data);
  return pdf->createNLL(dataHist)->getVal();
  pdf->fitTo(*data,RooFit::Extended(),RooFit::Strategy(0),RooFit::SumW2Error(kFALSE));
  RooFitResult* res = pdf->fitTo(*data,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2),RooFit::SumW2Error(kFALSE));
  std::cout << res->minNll() << std::endl;
  return res->minNll();
  */
}

TH1F* MakeSpinToy::getHistogram(RooAbsData* data, TString title, int rebin){
  TH1F* out = new TH1F(title,"",cosT->getBins()/rebin,cosT->getMin(),1);

  int i=0;
  while(data->get(i)){
    out->SetBinContent(i+1,data->weight());
    out->SetBinError(i+1,data->weightError());
    i++;
  }
  return out;
}

TTree* MakeSpinToy::makeForCombineTool(TString treeName, RooAbsData* hggData, RooAbsData* altData,RooAbsData* dataData){
  TTree * out = new TTree(treeName,"");
  float q;
  float gMf;
  int type;

  out->Branch("q",&q);
  out->Branch("NgenMinusNfit",&gMf);
  out->Branch("type",&type,"type/I");

  Long64_t iEntry=-1;
  type=1; //SM HIggs
  const RooArgSet *set;

  std::cout << "SM" <<std::endl;
  while( (set=hggData->get(++iEntry)) ){ 
    q = ((RooRealVar*)set->find("S"))->getVal();
    gMf = ((RooRealVar*)set->find("NgenMinusNfit"))->getVal();
    std::cout << q << std::endl;
    out->Fill();
  }
  iEntry=-1;
  type=-1;

  std::cout << "ALT" <<std::endl;
  while( (set=altData->get(++iEntry)) ){ 
    q = ((RooRealVar*)set->find("S"))->getVal();
    gMf = ((RooRealVar*)set->find("NgenMinusNfit"))->getVal();
    std::cout << q << std::endl;
    out->Fill();
  }
  if(dataData != 0){
  iEntry=-1;
  type=0;

  while( (set=dataData->get(++iEntry)) ){ 
    q = ((RooRealVar*)set->find("S"))->getVal();
    gMf = ((RooRealVar*)set->find("NgenMinusNfit"))->getVal();
    out->Fill();
  }
    
  }
  return out;
}

float MakeSpinToy::getExpEvents(float lumi, TString cat,TString mcType,RooWorkspace *toyws){
  //double tot    = ws->data(Form("%s_Combined",mcType.Data()))->sumEntries();
  //double thisN  = ws->data(Form("%s_%s",mcType.Data(),cat.Data()))->sumEntries();

  double f = ws->var(Form("Data_%s_FULLFIT_%s_fsig",mcType.Data(),cat.Data()))->getVal();
  
  if(toyws){
    toyws->import(*ws->var(Form("%s_EB_totalEvents",mcType.Data())));
    toyws->import(*ws->var(Form("%s_EE_totalEvents",mcType.Data())));
  }

  return f*lumi/12*607; //607 events in 12/fb @ 8 TeV
}


double* MakeSpinToy::run1(genType gen, int& N){
  TString mcType = mcLabels[gen];
  //we currently have 3 outputs
  N=7;
  double * out = new double[N];

  RooWorkspace *toyws=0;
  if(gen!=Data){
    toyws = new RooWorkspace(Form("toyws_%s",mcType.Data()),"toyws");
    generateToyWorkspace(toyws,gen);
    RooCategory* cosTCats = (RooCategory*)ws->obj("CosThetaBins");
    toyws->import(*cosTCats);

    MakeSpinFits fits("","");
    fits.setBkgFit(fitType);
    fits.setWorkspace(toyws);
    
    //float bins[6] = {0.0,0.2,0.4,0.6,0.8,1.0};
    //fits.setCosTBins(6,bins);

    fits.binDatasetCosT(*toyws->data("Data_Combined"),"Data");
    fits.binDatasetCosT(*toyws->data(mcLabels[1]+"_Combined"),mcLabels[1]);
    fits.binDatasetCosT(*toyws->data(mcLabels[2]+"_Combined"),mcLabels[2]);
    for( std::vector<TString>::const_iterator it = catLabels.begin();
	 it != catLabels.end(); it++){
      //toyws->import(*ws->data(Form("Hgg125_%s",it->Data())));
      //toyws->import(*ws->data(mcLabels[2]+"_"+*it));
      toyws->import(*ws->pdf(Form("%s_FIT_%s_cosTpdf",mcLabels[1].Data(),it->Data())));
      toyws->import(*ws->pdf(Form("%s_FIT_%s_cosTpdf",mcLabels[2].Data(),it->Data())));
      
      fits.MakeSignalFitForFit(*it,mcLabels[1]);
      fits.MakeSignalFitForFit(*it,mcLabels[2]);
      fits.MakeBackgroundOnlyFit(*it);          
      
      for(int iCos=0;iCos<cosTCats->numBins("");iCos++){
	cosTCats->setIndex(iCos);
	std::pair<float,float> cosTrange = MakeSpinFits::getCosTRangeFromCatName(cosTCats->getLabel());
	toyws->import(*ws->pdf( mcLabels[1]+"_FIT_"+*it+"_"+cosTCats->getLabel()) );
	toyws->import(*ws->pdf( mcLabels[2]+"_FIT_"+*it+"_"+cosTCats->getLabel()) );
	fits.MakeSignalFitForFit(*it+"_"+cosTCats->getLabel(),mcLabels[1]);
	fits.MakeSignalFitForFit(*it+"_"+cosTCats->getLabel(),mcLabels[2]);
	fits.MakeBackgroundOnlyFit(*it,cosTrange.first,cosTrange.second);
	}
    }
    toyws->import(*ws->pdf(Form("%s_FIT_cosTpdf",mcLabels[1].Data())));
    toyws->import(*ws->pdf(Form("%s_FIT_cosTpdf",mcLabels[2].Data())));

    //fits.MakeCombinedSignalTest(mcLabels[1]);
    //fits.MakeCombinedSignalTest(mcLabels[2]);
    //fits.MakeCombinedSignalTest(mcLabels[1],true); // make the binned signal tests
    //fits.MakeCombinedSignalTest(mcLabels[2],true);
    //fits.AddCombinedBkgOnlySWeight(mcLabels[1]);
    //fits.AddCombinedBkgOnlySWeight(mcLabels[2]);
    fits.MakeFullSBFit(mcLabels[1],false);          
    fits.MakeFullSBFit(mcLabels[1],true);          
    fits.MakeFullSBFit(mcLabels[2],false);          
    fits.MakeFullSBFit(mcLabels[2],true);          

    for( std::vector<TString>::const_iterator it = catLabels.begin();
	 it != catLabels.end(); it++){
      //fits.getSimpleBkgSubtraction(mcLabels[1],*it);
      //fits.getSimpleBkgSubtraction(mcLabels[2],*it);
    }
    
    //fits.Make2DCombinedSignalTest(mcLabels[2],mcLabels[2]);
    //fits.Make2DCombinedSignalTest(mcLabels[1],mcLabels[1]);

    fits.Make2DTemplateSignalTest(mcLabels[1]);
    fits.Make2DTemplateSignalTest(mcLabels[2]);

    //fits.getSimpleTotalBkgSubtraction(mcLabels[1]);
    //fits.getSimpleTotalBkgSubtraction(mcLabels[2]);

  }
  else{ //doing data
    toyws = ws;
    MakeSpinFits fits("","");
    fits.setWorkspace(toyws);
    fits.Make2DCombinedSignalTest(mcLabels[2],mcLabels[2]);
    fits.Make2DCombinedSignalTest(mcLabels[1],mcLabels[1]);
    fits.Make2DTemplateSignalTest(mcLabels[1]);
    fits.Make2DTemplateSignalTest(mcLabels[2]);

  }
  //toyws->Print();

  double hggll=0,altll=0;
  double shggll=0,saltll=0;
  double nFit=0,nGen=0;
  if(gen==Hgg125) nFit = toyws->var(Form("Data_%s_FULLSBFIT_Nsig",mcLabels[1].Data()))->getVal();
  else nFit = toyws->var(Form("Data_%s_FULLSBFIT_Nsig",mcLabels[2].Data()))->getVal();
  /*
  for( std::vector<TString>::const_iterator it = catLabels.begin();
       it != catLabels.end(); it++){
    //fits.MakeBackgroundFit("Hgg125",*it,125,2,false);

    RooDataHist *thisBkgHgg = (RooDataHist*)toyws->data( Form("Data_%s_%s_bkgSub_cosT",mcLabels[1].Data(),it->Data()));
    RooDataHist *thisBkgALT = (RooDataHist*)toyws->data( Form("Data_%s_%s_bkgSub_cosT",mcLabels[2].Data(),it->Data()));
    RooDataSet *thisSdataHgg = (RooDataSet*)toyws->data( Form("Data_%s_%s_sigWeight",it->Data(),mcLabels[1].Data()));
    RooDataSet *thisSdataALT = (RooDataSet*)toyws->data( Form("Data_%s_%s_sigWeight",it->Data(),mcLabels[2].Data()));
    RooDataHist *thisShistHgg = new RooDataHist("tmpHgg","",*cosT,*thisSdataHgg);
    RooDataHist *thisShistALT = new RooDataHist("tmpALT","",*cosT,*thisSdataALT);

    hggPdf = (RooHistPdf*)ws->pdf(Form("%s_FIT_%s_cosTpdf",mcLabels[1].Data(),it->Data()));
    altPdf = (RooHistPdf*)ws->pdf(Form("%s_FIT_%s_cosTpdf",mcLabels[2].Data(),it->Data()));


    if(gen!=Data) nGen+=toyws->var( Form("N_gen_sig_%s",it->Data()) )->getVal();
    else nGen+=thisSdataHgg->sumEntries();

    if(saveWorkspaces && gen!=Data){
      toyws->import(*hggPdf);
      toyws->import(*altPdf);
    }
    
    std::cout << "hgg" <<std::endl;
    hggll += computeLL(hggPdf,thisBkgHgg,cosT,2);
    std::cout << "alt" <<std::endl;
    //altll += computeLL(altPdf,thisBkgALT,cosT,3);
    altll += computeLL(altPdf,thisBkgHgg,cosT,2);
    
    std::cout << "hgg splot" <<std::endl;
    shggll += computeLL(hggPdf,thisShistHgg,cosT,2);
    std::cout << "alt splot" <<std::endl;
    //saltll += computeLL(altPdf,thisShistALT,cosT,3);
    saltll += computeLL(altPdf,thisShistHgg,cosT,2);

    std::cout << "Bkg-sub:" << std::endl;
    std::cout << "Hgg: " << hggll << std::endl;
    std::cout << "ALT: " << altll << std::endl;

    std::cout << "splot:" << std::endl;
    std::cout << "Hgg: " << shggll << std::endl;
    std::cout << "ALT: " << saltll << std::endl;

    delete thisShistHgg;
    delete thisShistALT;
  }
  */
  /*
  RooDataHist *thisBkgHgg = (RooDataHist*)toyws->data(Form("Data_%s_Combined_bkgSub_cosT",mcLabels[1].Data()));
  RooDataHist *thisBkgALT = (RooDataHist*)toyws->data(Form("Data_%s_Combined_bkgSub_cosT",mcLabels[2].Data()));
  hggPdf = (RooHistPdf*)ws->pdf(mcLabels[1]+"_FIT_cosTpdf");
  altPdf = (RooHistPdf*)ws->pdf(mcLabels[2]+"_FIT_cosTpdf");

  double totHggLL = computeLL(hggPdf,thisBkgHgg,cosT);
  //double totALTLL = computeLL(altPdf,thisBkgALT,cosT);
  double totALTLL = computeLL(altPdf,thisBkgHgg,cosT);

  //use 2D fit to test compatibility
  RooRealVar *yield1D   = toyws->var(Form("Data_%s_FULLFIT_Nsig",mcLabels[1].Data()));
  RooRealVar *yieldHgg  = toyws->var(Form("Data_m_%s_c_%s_FULL2DFIT_Nsig",mcLabels[1].Data(),mcLabels[1].Data())); // yield fitting hgg mass and hgg cosT distributions
  RooRealVar *yieldALT  = toyws->var(Form("Data_m_%s_c_%s_FULL2DFIT_Nsig",mcLabels[2].Data(),mcLabels[2].Data())); // yield fitting alt mass and alt cosT distributions
  
  double errHgg = TMath::Sqrt(pow(yield1D->getError(),2) + pow(yieldHgg->getError(),2));
  double errALT = TMath::Sqrt(pow(yield1D->getError(),2) + pow(yieldALT->getError(),2));

  double hgg2Dll  = TMath::Log(TMath::Prob( pow(yield1D->getVal()-yieldHgg->getVal(),2)/pow(errHgg,2),1));
  double alt2Dll  = TMath::Log(TMath::Prob( pow(yield1D->getVal()-yieldALT->getVal(),2)/pow(errALT,2),1));

  RooFitResult *fitHgg = (RooFitResult*)toyws->obj(Form("Data_m_%s_c_%s_FULL2DFIT_fitResult",mcLabels[1].Data(),mcLabels[1].Data())); 
  RooFitResult *fitALT = (RooFitResult*)toyws->obj(Form("Data_m_%s_c_%s_FULL2DFIT_fitResult",mcLabels[2].Data(),mcLabels[2].Data()));
  */
  RooFitResult *fitTempHgg = (RooFitResult*)toyws->obj(Form("Data_%s_TEMPLATE2DFIT_fitResult",mcLabels[1].Data()));
  RooFitResult *fitTempALT = (RooFitResult*)toyws->obj(Form("Data_%s_TEMPLATE2DFIT_fitResult",mcLabels[2].Data()));

  RooFitResult *fitCostHgg = (RooFitResult*)toyws->obj(Form("Data_%s_FULLSBCOSTFIT_fitResult",mcLabels[1].Data()));
  RooFitResult *fitCostALT = (RooFitResult*)toyws->obj(Form("Data_%s_FULLSBCOSTFIT_fitResult",mcLabels[2].Data()));

  
  /*
  std::cout << "FIT NLL: " << endl
	    << "Hgg: " << fitHgg->minNll() << endl
	    << "ALT: " << fitALT->minNll() << endl;
  double hgg2DFITNll = fitHgg->minNll();
  double alt2DFITNll = fitALT->minNll();
  */
  double hgg2DTEMPFITNll = fitTempHgg->minNll();
  double alt2DTEMPFITNll = fitTempALT->minNll();
  
  double hggCostFITNll = fitCostHgg->minNll();
  double altCostFITNll = fitCostALT->minNll();
  std::cout << "BIN FIT NLL: " << endl
	    << "Hgg: " << fitCostHgg->minNll() << endl
	    << "ALT: " << fitCostALT->minNll() << endl;
  

  GenMinusFit->setVal( nGen-nFit);
  std::cout << "Gen - Fit: " << nGen - nFit <<std::endl;
  if(saveWorkspaces) toyWSs.AddLast(toyws);
  else delete toyws;

  //setup outputs:
  /*
  out[0] = 2*(altll-hggll);
  out[1] = 2*(totALTLL-totHggLL);
  out[2] = 2*(saltll-shggll);
  out[3] = 2*(alt2Dll-hgg2Dll);
  out[4] = -2*(alt2DFITNll-hgg2DFITNll);
  */
  out[0]=out[1]=out[2]=out[3]=out[4]=0;
  out[5] = -2*(alt2DTEMPFITNll-hgg2DTEMPFITNll);
  out[6] = -2*(altCostFITNll-hggCostFITNll);

  return out;
}

void MakeSpinToy::runN(int N){
  RooMsgService::instance().setSilentMode(true);
  int Nout;
  for(int i=0;i<N;i++){
    double * hggOut = run1(Hgg125,Nout);

    S->setVal(hggOut[0]);
    S_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[1]);
    S_tot_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[2]);
    S_splot_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[3]);
    S_2D_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[4]);
    S_2DFIT_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[5]);
    S_2DTEMPFIT_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[6]);
    S_2DBINFIT_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    double * altOut = run1(ALT125,Nout);

    S->setVal(altOut[0]);
    S_TruthALT->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(altOut[1]);
    S_tot_TruthALT->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(altOut[2]);
    S_splot_TruthALT->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(altOut[3]);
    S_2D_TruthALT->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(altOut[4]);
    S_2DFIT_TruthALT->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(altOut[5]);
    S_2DTEMPFIT_TruthALT->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(altOut[6]);
    S_2DBINFIT_TruthALT->add(RooArgSet(*S,*GenMinusFit));   

    delete hggOut;
    delete altOut;
  }

  if(doData){
    double *dataOut = run1(Data,Nout);
    S->setVal(dataOut[0]);
    S_TruthData->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(dataOut[1]);
    S_tot_TruthData->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(dataOut[2]);
    S_splot_TruthData->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(dataOut[3]);
    S_2D_TruthData->add(RooArgSet(*S,*GenMinusFit));   
    
    S->setVal(dataOut[4]);
    S_2DFIT_TruthData->add(RooArgSet(*S,*GenMinusFit));   
    
    S->setVal(dataOut[5]);
    S_2DTEMPFIT_TruthData->add(RooArgSet(*S,*GenMinusFit));   
    
    S->setVal(dataOut[6]);
    S_2DBINFIT_TruthData->add(RooArgSet(*S,*GenMinusFit));   
    
  }
}

void MakeSpinToy::save(TString outputFile){
  TFile *f = new TFile(outputFile,"RECREATE");
  makeForCombineTool("q",S_TruthHgg,S_TruthALT,S_TruthData)->Write();
  makeForCombineTool("qtot",S_tot_TruthHgg,S_tot_TruthALT,S_tot_TruthData)->Write();
  makeForCombineTool("qsplot",S_splot_TruthHgg,S_splot_TruthALT,S_splot_TruthData)->Write();
  makeForCombineTool("q2D",S_2D_TruthHgg,S_2D_TruthALT,S_2D_TruthData)->Write();
  makeForCombineTool("q2DFIT",S_2DFIT_TruthHgg,S_2DFIT_TruthALT,S_2DFIT_TruthData)->Write();
  makeForCombineTool("q2DTEMPFIT",S_2DTEMPFIT_TruthHgg,S_2DTEMPFIT_TruthALT,S_2DTEMPFIT_TruthData)->Write();
  makeForCombineTool("q2DBINFIT",S_2DBINFIT_TruthHgg,S_2DBINFIT_TruthALT,S_2DBINFIT_TruthData)->Write();
  if(toyWSs.GetEntries()>0) toyWSs.Write();
  f->Close();
}
