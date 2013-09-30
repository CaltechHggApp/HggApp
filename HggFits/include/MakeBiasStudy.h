#ifndef MakeBiasStudy_h
#define MakeBiasStudy_h
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
#include "RooArgSet.h"
#include "RooLinkedListIter.h"
#include "RooAddPdf.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"

#include "MakeSpinFits.h"
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

#define NUM_CPU 1

class MakeBiasStudy {
public:
    MakeBiasStudy(const TString& inputFileName, const TString& outputFileName); //!< constructor requires file path info
    virtual ~MakeBiasStudy();
    /*!
      make a fake signal model from a previous one scaling the widths appropriately and setting the mean
    */

    void biasStudy(MakeSpinFits::BkgFitType truthType, const int N); //!< make the bias study
    void biasStudy(MakeSpinFits::BkgFitType truthType,const int N, const TString& cat); //!< make the bias study for a given category

    std::tuple<double,double,double> getNFit(RooAbsData& toyData,MakeSpinFits::BkgFitType fitType, float sigma_eff,RooAbsPdf& signalModel);

    void run();
    void print();

    void setNToys(int n){NToys=n;}
    
    typedef std::map<TString,std::map<MakeSpinFits::BkgFitType,std::vector<double> > > biasList;
private:
    RooWorkspace *ws;

    TFile* inputFile=0;
    TFile* outputFile=0;
    std::vector<TString> catLabels;

    std::vector<MakeSpinFits::BkgFitType> testFitTypes  = {MakeSpinFits::kPow,MakeSpinFits::kDoubleExp,MakeSpinFits::kModifiedExp,MakeSpinFits::kPoly};
    std::vector<MakeSpinFits::BkgFitType> truthFitTypes = {MakeSpinFits::kPow,MakeSpinFits::kDoubleExp,MakeSpinFits::kModifiedExp};
    //std::vector<MakeSpinFits::BkgFitType> testFitTypes  = {MakeSpinFits::kPow,MakeSpinFits::kDoubleExp};
    //std::vector<MakeSpinFits::BkgFitType> truthFitTypes = {MakeSpinFits::kPow};
    biasList biasSigError;
    biasList biasBkgError;
    float mh=125;

    int NToys = 100;

    std::map<MakeSpinFits::BkgFitType,TString> fitNameMap = {
        {MakeSpinFits::kPow, "Single Power Law"}, {MakeSpinFits::kDoublePow, "Double Power Law"}, {MakeSpinFits::kSingleExp,"Single Exponential"},
        {MakeSpinFits::kDoubleExp, "Double Exponential"},{MakeSpinFits::kTripleExp,"Triple Exponential"},{MakeSpinFits::kModifiedExp,"Modified Exponential"},
        {MakeSpinFits::kPoly,"5th Order Polynomial"}                                                                                                          
    };

    void printFormatted(const biasList& list);
};

#endif
