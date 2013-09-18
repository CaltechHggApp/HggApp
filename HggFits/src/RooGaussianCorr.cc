
#include "Riostream.h"

#include "include/RooGaussianCorr.hh"
#include "TMath.h"

ClassImp(RooGaussianCorr)


RooGaussianCorr::RooGaussianCorr(const char* name, const char *title, RooArgList &variables, RooArgList &means, TMatrixDSym *covarianceMatrix) :
  RooAbsPdf(name,title),
  __variables("vars","variable list",this),
  __means("means","mean list",this),
  __covMatrix(covarianceMatrix),
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
  __covMatrix(other.__covMatrix),
  __invCovMatrix(0),
  determinant(other.determinant)
{
  invert();
}

void RooGaussianCorr::invert() { //invert the covariance matrix
  if(__invCovMatrix) delete __invCovMatrix;
  __invCovMatrix = new TMatrixDSym(*__covMatrix); //create a new object which we own
  __invCovMatrix->Invert(); //invert in place
}

double RooGaussianCorr::evaluate() const {
  assert(__invCovMatrix != 0);

  int k = __invCovMatrix->GetNcols(); // the asserts ensure that this is the same dimension as the others
  double norm = 1./TMath::Sqrt(TMath::Power(TMath::TwoPi(),k)*determinant); //normalization factor

  TVectorD varMinusMean(k);
  for(int i=0;i<k;i++) { // build the vector in place
    RooAbsReal *var  = (RooAbsReal*)__variables.at(i);
    RooAbsReal *mean = (RooAbsReal*)__means.at(i);
    varMinusMean(i) = var->getVal() - mean->getVal();
  }
  TVectorD varMinusMean2 = varMinusMean; // do the multiplication (v.I.v)
  varMinusMean *= *__invCovMatrix;       // we need some copying, but this should minimize it
  double expfactor = varMinusMean * varMinusMean2;
  return norm*TMath::Exp(-0.5*expfactor);
}
