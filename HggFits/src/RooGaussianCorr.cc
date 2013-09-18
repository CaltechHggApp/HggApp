
#include "Riostream.h"

#include "include/RooGaussianCorr.hh"
#include "TMath.h"

ClassImp(RooGaussianCorr)


RooGaussianCorr::RooGaussianCorr(const char* name, const char *title, RooArgList &variables, RooArgList &means, const TMatrixDSym *covarianceMatrix) :
  RooAbsPdf(name,title),
  __variables("vars","variable list",this),
  __means("means","mean list",this),
  __covMatrix(new TMatrixDSym(*covarianceMatrix)),
  __invCovMatrix(0)
{
  //build the fing RooListProxies
  TIterator* _varIter  = variables.createIterator();
  TIterator* _meanIter = means.createIterator();

  RooAbsReal *var;
  RooAbsReal *mean;
  while( (var=(RooAbsReal*)_varIter->Next()) ) {
    mean = (RooAbsReal*)_meanIter->Next();
    if(!mean) {
      std::cout << "number of vars != number of means" << std::endl;
      assert(0);
    }
    __variables.add(*var);
    __means.add(*mean);
  }
  _varIter->Next();
  if(var) {
    std::cout << "number of vars != number of means" << std::endl;
    assert(0);
  }

  determinant = __covMatrix->Determinant();
  assert(__covMatrix->GetNcols() == __variables.getSize());
  assert(__covMatrix->GetNcols() == __means.getSize());
  invert();
}

RooGaussianCorr::RooGaussianCorr(const RooGaussianCorr & other, const char* name) :
  RooAbsPdf(other,name),
  __variables("vars",this,other.__variables),
  __means("means",this,other.__means),
  __covMatrix(new TMatrixDSym(*(other.__covMatrix.get()))),
  __invCovMatrix(0),
  determinant(other.determinant)
{
  invert();
}

void RooGaussianCorr::invert() { //invert the covariance matrix
    if(isInverted) return;
    isInverted=true;
    if(__invCovMatrix) delete __invCovMatrix;
    assert(__covMatrix!=0);
    __invCovMatrix = (TMatrixDSym*)__covMatrix->Clone("invertedCovarianceMatrix"); //new TMatrixDSym(*__covMatrix); //create a new object which we own
    __invCovMatrix->Invert(); //invert in place

    dim = __invCovMatrix->GetNcols(); // the asserts ensure that this is the same dimension as the others
}
#include <iostream>
double RooGaussianCorr::evaluate() const {
  assert(__invCovMatrix != 0);
  assert(isInverted);
  assert(dim>0);
  assert(determinant!=0);  //evaluate will be deep in the guts of the fitter, so use assert to stop immediately
  
  //double norm = 1./TMath::Sqrt(TMath::Power(TMath::TwoPi(),dim)*determinant); //normalization factor
  double norm = 1; //getVal() takes care of the normalization
  assert(norm!=0);
  TVectorD varMinusMean(dim);
  for(int i=0;i<dim;i++) { // build the vector in place
    RooAbsReal *var  = (RooAbsReal*)__variables.at(i);
    RooAbsReal *mean = (RooAbsReal*)__means.at(i);
    varMinusMean(i) = var->getVal() - mean->getVal();
  }
  TVectorD varMinusMean2 = varMinusMean; // do the multiplication (v.I.v)
  varMinusMean *= *__invCovMatrix;       // we need some copying, but this should minimize it
  double expfactor = varMinusMean * varMinusMean2;
  double ret_val = norm*TMath::Exp(-0.5*expfactor);

  /*
  if(ret_val < 1e-20){
      std::cout << ret_val  << std::endl;
      varMinusMean2.Print();
      __invCovMatrix->Print();
      varMinusMean.Print();
  }
  */
  return norm*TMath::Exp(-0.5*expfactor);

}
