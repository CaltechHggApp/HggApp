#include "MakeBiasStudy.h"

#include "assert.h"



MakeBiasStudy::MakeBiasStudy(const TString& inputFileName, const TString& outputFileName)
{
  if(inputFileName != ""){
    //opens the input file and gets the input workspace
      inputFile = new TFile(inputFileName);
      ws = ((RooWorkspace*)inputFile->Get("cms_hgg_spin_workspace"));
      //ws->Print();
      //extract MC labels from the input workspace
      //MakeSpinFits::getLabels("labels",&mcLabel,ws);
      //extract category labels from the input workspace
      MakeSpinFits::getLabels("evtcat",&catLabels,ws);
  }
  if(outputFileName != ""){
      //opens the output file
      outputFile = new TFile(outputFileName,"RECREATE");
    outputFile->cd();
    ws->Write(ws->GetName(),TObject::kWriteDelete);
  }

}

MakeBiasStudy::~MakeBiasStudy() {
}

void MakeBiasStudy::biasStudy(MakeSpinFits::BkgFitType truthType,const int N) {
    for(auto cat: catLabels) {
        biasStudy(truthType,N,cat);        
    }
}

void MakeBiasStudy::biasStudy(MakeSpinFits::BkgFitType truthType, const int N, const TString& cat) {
    RooRealVar* mass = ws->var("mass");

    //get real data and fit a pdf of type truthType to it
    RooAbsData* realData = ws->data("Data_Combined")->reduce( Form("evtcat==evtcat::%s",cat.Data()) );
    RooAbsPdf *truthPdf  = MakeSpinFits::getBackgroundPdf("Truth","BKGFIT",cat,truthType,mass);
    RooRealVar nBkgTruth("TruthNBkg","",0,1e9);
    RooExtendPdf truthExtendedPdf("truthExtendedPdf","",*truthPdf,nBkgTruth);
    truthExtendedPdf.fitTo(*realData,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    truthExtendedPdf.fitTo(*realData,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

    //get the signal model for this category and compute its sigma effective
    RooAbsPdf* truthSignal = ws->pdf(  Form("jhu0plus125_FIT_%s",cat.Data())  );
    assert(truthSignal!=0);
    double sigEff = MakeSpinFits::computeSigEff(truthSignal,mh,mass);
    
    RooRealVar range("range","",mh-sigEff,mh+sigEff);
    RooRealVar all("all","",ws->var("mass")->getMin(),ws->var("mass")->getMax());

    double BkgInRange  = truthExtendedPdf.createIntegral(range)->getVal()/truthExtendedPdf.createIntegral(all)->getVal()*nBkgTruth.getVal();
    double BkgInRangeE = truthExtendedPdf.createIntegral(range)->getVal()/truthExtendedPdf.createIntegral(all)->getVal()*nBkgTruth.getError();

    std::cout << "TTTT" << BkgInRange << std::endl;


    //store results for the signal and background error 
    std::vector<std::vector<double>> biasesSigErr(testFitTypes.size());
    std::vector<std::vector<double>> biasesBkgErr(testFitTypes.size());

    /*
      for(int i=0;i<testFitTypes.size();i++) {
        biasesSigErr.push_back(std::vector<double>>(new std::vector<double>));
        biasesBkgErr.push_back(new std::vector<double>));
    }
    */

    //throw toys
    for(int i=0;i<N;i++) {
      //throw toy background + signal data
        RooDataSet * toyData = truthExtendedPdf.generate(RooArgSet(*mass),int(nBkgTruth.getVal()),RooFit::Extended(true));
        toyData->append(*truthSignal->generate(RooArgSet(*mass),
                                               int(ws->data("jhu0plus125_Combined")->reduce(Form("evtcat==evtcat::%s",cat.Data()))->sumEntries()),
                                               RooFit::Extended(true)));
        
        //ws->import(*toyData,RooFit::Rename(Form("toy_%s",cat.Data())));
        auto fitIt = testFitTypes.begin();
        auto biasSEIt = biasesSigErr.begin();
        auto biasBEIt = biasesBkgErr.begin();
	//for each testFitType, compute that uncertainty on the background and signal yield
        for(; fitIt!=testFitTypes.end();fitIt++,biasSEIt++, biasBEIt++) {            
            std::tuple<double,double,double> results = getNFit(*toyData,*fitIt,sigEff,*truthSignal);
            biasSEIt->push_back( (std::get<0>(results)-BkgInRange)/std::get<2>(results) );
            biasBEIt->push_back( (std::get<0>(results)-BkgInRange)/std::get<1>(results) );
        }
    }
    auto fitIt = testFitTypes.begin();
    auto biasSEIt = biasesSigErr.begin();
    auto biasBEIt = biasesBkgErr.begin();
    for(; fitIt!=testFitTypes.end();fitIt++,biasSEIt++, biasBEIt++) {            
        std::cout << *fitIt << std::endl;
        std::cout << biasesSigErr.size() << std::endl;
        std::cout << biasSEIt->size() << std::endl;
        std::sort( biasSEIt->begin(), biasSEIt->end() );
        std::sort( biasBEIt->begin(), biasBEIt->end() );       

	//get the median error
        biasSigError[cat][truthType].push_back( biasSEIt->at( N/2 ));
        biasBkgError[cat][truthType].push_back( biasBEIt->at( N/2 ));
        std::cout << "QQQQQQQQQQQ BIAS SIG ERROR Truth: " << fitNameMap[truthType] << "   Fit: " << fitNameMap[*fitIt] << std::endl;
        for(auto bias : *biasSEIt) {
            std::cout << bias << " ";
        }
        std::cout << std::endl;
        std::cout << "QQQQQQQQQQQ BIAS BKG ERROR Truth: " << fitNameMap[truthType] << "   Fit: " << fitNameMap[*fitIt] << std::endl;
        for(auto bias : *biasBEIt) {
            std::cout << bias << " ";
        }
        std::cout << std::endl;
    }
    //ws->import(*truthPdf);
    delete truthPdf;
}

std::tuple<double,double,double> MakeBiasStudy::getNFit(RooAbsData& toyData, MakeSpinFits::BkgFitType fitType, float sigma_eff,RooAbsPdf& signalModel) {
    RooRealVar* mass = ws->var("mass");


    //build a pdf of type fitType and fit it to the toyData
    RooAbsPdf *fitPdf = MakeSpinFits::getBackgroundPdf("FIT","BKGFIT","f",fitType,mass);
    RooRealVar nBkgFit("FitNBkg","",0,1e9);
    RooExtendPdf fitExtendedPdf("fitExtendedPdf","",*fitPdf,nBkgFit);
    fitExtendedPdf.fitTo(toyData,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    fitExtendedPdf.fitTo(toyData,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

    //compute the amount of background in the +- 1 sigma_eff range
    RooRealVar range("range","",mh-sigma_eff,mh+sigma_eff);
    RooRealVar all("all","",ws->var("mass")->getMin(),ws->var("mass")->getMax());

    double BkgInRange  = fitExtendedPdf.createIntegral(range)->getVal()/fitExtendedPdf.createIntegral(all)->getVal()*nBkgFit.getVal();
    double BkgInRangeE = fitExtendedPdf.createIntegral(range)->getVal()/fitExtendedPdf.createIntegral(all)->getVal()*nBkgFit.getError();

    std::cout << "NNNN" << BkgInRange << std::endl;
    
    //fix all the background parameters
    RooArgSet *bkg_pars = fitExtendedPdf.getVariables();
    RooFIter iter = bkg_pars->fwdIterator();
    RooAbsArg* a;
    while( (a = iter.next()) ){
        if(string(a->GetName()).compare("mass")==0) continue;
        static_cast<RooRealVar*>(a)->setConstant(kTRUE);
    }

    //build a S+B pdf with the fixed background model just computed and the signalModel provided
    RooRealVar nBkgSB("NBkgSB","",0,1e9);
    RooRealVar nSigSB("NSigSB","",0,1e5);
    
    RooAddPdf sumPdf("sumPdf","",RooArgSet(fitExtendedPdf,signalModel),RooArgSet(nBkgSB,nSigSB));
    sumPdf.fitTo(toyData,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    sumPdf.fitTo(toyData,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

    delete fitPdf;
    delete bkg_pars;
    return std::tuple<double,double,double>(BkgInRange,BkgInRangeE,nSigSB.getError());
}

            
void MakeBiasStudy::run() {

    for(auto truthIt : truthFitTypes) {
        biasStudy(truthIt,NToys);
    }
    print();
    outputFile->cd();
    //ws->Write();
    outputFile->Close();
}

void MakeBiasStudy::print() {
    std::cout << "\n\n*****************SIGNAL ERROR********************\n" << std::endl;
    printFormatted(biasSigError);
    std::cout << "\n\n*****************BKG    ERROR********************\n" << std::endl;
    printFormatted(biasBkgError);
}
void MakeBiasStudy::printFormatted(const biasList& list) {
    const char *form = "%40s";
    
    for(auto catIt : catLabels) {
        std::cout << catIt << std::endl;
        std::cout << Form(form,"");
        for(auto fitIt: testFitTypes) {
            std::cout << Form(form,fitNameMap[fitIt].Data());
        }
        std::cout << std::endl;
        for(auto fitIt: truthFitTypes) {
            std::cout << Form(form,fitNameMap[fitIt].Data());
            for(auto biasIt : list.at(catIt).at(fitIt) ){
                std::cout << Form("%40.2f",biasIt);
            }
            std::cout << std::endl;
        }
    }
}

