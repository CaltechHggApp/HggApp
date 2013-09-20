#ifndef MakeSimpleHypothesisTest_h
#define MakeSimpleHypothesisTest_h
//! class to do simple hypothesis testing over a mass range

/*!

This class takes an input workspace with fits and does simple signal testing over the mass range of the mass variables

Author: Alex Mott (Caltech)
Date: Sep 2014
*/

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"

#include "MakeSpinFits.h"

class MakeSimpleHypothesisTest {
public:
    MakeSimpleHypothesisTest(const TString& inputFileName, const TString& outputFileName); //!< constructor requires file path info

    /*!
      make a fake signal model from a previous one scaling the widths appropriately and setting the mean
    */
    void MakeFakeSignalModels(const TString& baseMCName, const TString& tag,float scaleFactor,float meanVal);
    void MakeAllFakeSignalModels();

    void MakeHypothesisTest();
    void run();
private:
    RooWorkspace *ws;

    TFile* inputFile=0;
    TFile* outputFile=0;
    std::vector<TString> catLabels;

    int mh_low=100;
    int mh_high=4900;
    int mh_step=1;
    
    
};

#endif
