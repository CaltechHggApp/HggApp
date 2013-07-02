#ifndef MakeSpinSPlot_h
#define MakeSpinSPlot_h

/*!
  A Generic Class for computing SWeights for a data set.
  The yields must be fitted externally and the fitted yields passed to the class

  The output will be a dataset that can be merged with the input datasets
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
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"

#include "TMatrixT.h"
#include "TMath.h"

class MakeSpinSPlot{
public:

  MakeSpinSPlot(RooAbsData *data); //!< Constructor, dataset is passed at this point
  
  void addSpecies(TString name, RooAbsPdf* pdf, double expYield); //!< add a species to have an sweight computed.  
  /*!
    \param name The name to use for the species (e.g. 'signal' and 'background') The SWeight will be named <name>_sw
    \param pdf  The pdf associated with this species
    \\param expYield the fitted yield of this species
  */
  void addVariable(RooRealVar* var){ __variables.push_back(var);} //!< add a variable
  /*!
    Add a variable whose value will be updated for every event.  All varaibles used in any species pdf must be given here.
  */
  
  void calculate(); //!< run the sweight calculation

  RooDataSet* getSWeightDataSet(){ return __sWeightDataSet; } //!< return the SWeight dataset

  RooArgSet* getSWeightVars(){return __sWeightVars;} //!< get the sweight varaibles
protected:
  void computeCovMatrix(); //!< internal method to compute the covariance matrix

  double computeDenom(); //!< compute the denominator for the covariance matrix and the SWeight

  void computeSWeight(); //!< internal computer
  RooDataSet* __sWeightDataSet;
  RooArgSet* __sWeightVars;

  int __nSpec;
  std::vector<RooAbsPdf*> __pdfs;
  std::vector<double> __expectedYields;
  std::vector<TString> __speciesNames;

  std::vector<RooRealVar*> __variables;

  RooAbsData *__dataSet;
  
  
  TMatrixD *__covMatrix;
};

#endif
