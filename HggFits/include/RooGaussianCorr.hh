//
// A multidimensional gaussian pdf that implements correlations amongst the input variables
// use-case is for a multi-dimensional constraint of fitting parameters
//

#include "RooAbsArg.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "TObject.h"

#include "RooListProxy.h"
#include "RooRealProxy.h"

#include "TMatrixDSym.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TIterator.h"

#include "assert.h"
#include <iostream>


class RooGaussianCorr : public RooAbsPdf {
public:
  RooGaussianCorr(const char* name, const char *title,RooArgList &variables, RooArgList &means, TMatrixDSym &covarianceMatrix);
  RooGaussianCorr(const RooGaussianCorr & other,const char* name=0);
  ~RooGaussianCorr() { delete __invCovMatrix; }
protected:
  RooListProxy __variables;
  RooListProxy __means;
  TMatrixDSym &__covMatrix;
  TMatrixDSym *__invCovMatrix;

  double determinant; // compute this once and save the cached result
  double evaluate() const;
private:
  void invert();
  ClassDef(RooGaussianCorr,1)
};
