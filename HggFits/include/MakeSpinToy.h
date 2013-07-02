#ifndef MakeSpinToy_h
#define MakeSpinToy_h
//! class to run spin toys on the fitted Hgg workspace

/*!

This class takes a workspace that has had all fits run on it and throws toys to extract
the spin separation. between two MC samples.

Author: Alex Mott (Caltech)
Date: Feb 2013
*/


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
#include <TObjArray.h>
#include <TH1F.h>
#include <math.h>

#include "MakeSpinFits.h"
#include "MakeSpinPlots.h"
#include "MakeSpinWorkspace.h"

class MakeSpinToy{
public:
  MakeSpinToy(TString fileName,TString wsName="cms_hgg_spin_workspace"); //!< default constructor takes the path to the workspace file
  void runMCStudy(int Ntoys,float lumi,TString cat);

  void setTargetLumi(float l){targetLumi = l;}
  void setNominalLumi(float l){nominalLumi = l;}

  enum genType{Data,Hgg125,ALT125};
  
  float getExpEvents(float lumi, TString cat,TString mcType,RooWorkspace* toyws=0); //!< get the expected number of higgs events for the given lumi and category

  /*!
    this function computes the binned p-value of the pdf to the given dataset and converts that into a log-likelihood

    \param pdf  a pointer to the pdf to test
    \param data a pointer to the data to test
    \param var  a pointer to the variable of the pdf
    \param rebin specify the number of bins to merge for the p-value computation.  The baseline number of bins is the number in var->getBins("")
  */
  double computeLL(RooAbsPdf* pdf, RooAbsData* data,RooRealVar* var, int rebin=1); //!< compute the LL of the PDF to the data

  /*! 
    this function runs a single toy workspace.  One specifies the truth model for the toy with the input genType.
    The function returns the separation computed in a number of different ways (N will be set to this number).
    the return value is a pointer to an array of size N containing the computed q values.
  */
  double* run1(genType gen, int& N); //!< run a single toy 

  void runN(int N); //!< run N toys with both SM and alternate truth models

  void save(TString outputFile); //!< save all the output trees

  TH1F* getHistogram(RooAbsData* data, TString title,int rebin=1); //!< make a TH1F out of the given data


  //these are the necessary workspaces, variables, datasets and pdfs for the calculations
  RooWorkspace *ws;
  RooRealVar* cosT, *mass;
  RooRealVar *S,*GenMinusFit;
  RooDataSet *S_TruthHgg, *S_TruthALT, *S_TruthData;
  RooDataSet *S_tot_TruthHgg, *S_tot_TruthALT, *S_tot_TruthData;
  RooDataSet *S_splot_TruthHgg, *S_splot_TruthALT, *S_splot_TruthData;
  RooDataSet *S_2D_TruthHgg, *S_2D_TruthALT, *S_2D_TruthData;
  RooDataSet *S_2DFIT_TruthHgg, *S_2DFIT_TruthALT, *S_2DFIT_TruthData;
  RooDataSet *S_2DTEMPFIT_TruthHgg, *S_2DTEMPFIT_TruthALT, *S_2DTEMPFIT_TruthData;
  RooDataSet *S_2DBINFIT_TruthHgg, *S_2DBINFIT_TruthALT, *S_2DBINFIT_TruthData;
  RooHistPdf *altPdf,*hggPdf,*bkgPdf;

  /*!
    this function takes in a dataset containing the separation for a set of toys and converts it to a TTree in a 
    format that the combine tool can read for final separation testing and plotting
  */
  TTree* makeForCombineTool(TString treeName, RooAbsData* hggData, RooAbsData* altData,RooAbsData* dataData=0); //!< convert the output to a format readable by the combine tool for plot making

  /*!
    The function tells the save function whether to save the toy workspaces.  By default the workspaces are deleted after the separation is extracted.  
    WARNING: the workspaces can be rather larger (tens of MB each), so this option should only be used for debugging.
  */
  void setSaveWorkspaces(bool b){saveWorkspaces = b;} //!< specify whether to save the toy workspaces

  void setDoData(bool b); //!< specify that we should also extract the q value from the real data (not just toys).

  void setMCStandard(TString s){mcLabels[1]=s;} //!< specify the name of the SM hypothesis (default: Hgg125)
  void setMCComparison(TString s){mcLabels[2]=s;} //!< specify the name of the alternate hypothesis

  void setBkgFitType(MakeSpinFits::BkgFitType f){fitType=f;} //!< set the type of background fit for the toys

  TObjArray toyWSs;
protected:
  float targetLumi;
  float nominalLumi;
  void generateToyWorkspace(RooWorkspace* toyws, const char* cat,genType gen,float nSigTot);
  void generateToyWorkspace(RooWorkspace* toyws,genType gen);
  
  bool useR9;
  bool saveWorkspaces;
  bool doData;
  std::vector<TString> catLabels;
  std::vector<TString> catCosTLabels;

  MakeSpinFits::BkgFitType fitType;
  
  TString mcLabels[3];
};

#endif
