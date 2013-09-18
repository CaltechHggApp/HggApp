#ifndef RooGaussianCorr_hh
#define RooGaussianCorr_hh

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
#include "TMemberInspector.h"
#include "TBuffer.h"

#include "assert.h"
#include <iostream>
#include <memory>

class RooGaussianCorr : public RooAbsPdf {
public:
    RooGaussianCorr() {};    
    RooGaussianCorr(const char* name, const char *title,RooArgList &variables, RooArgList &means, const TMatrixDSym *covarianceMatrix);
    RooGaussianCorr(const TString name, const char *title,RooArgList &variables, RooArgList &means, const TMatrixDSym *covarianceMatrix) :
        RooGaussianCorr(name.Data(),title,variables,means,covarianceMatrix) {}
    RooGaussianCorr(const RooGaussianCorr & other,const char* name=0);
    virtual TObject* clone(const char* newname) const { return new RooGaussianCorr(*this,newname); }
    inline virtual ~RooGaussianCorr() { delete __invCovMatrix; }
    
  //virtual void ShowMembers(TMemberInspector& insp, char* parent){}
  //virtual void Streamer(TBuffer& b) {}
protected:
    RooListProxy __variables;
    RooListProxy __means;
    std::unique_ptr<TMatrixDSym> __covMatrix; //copy the covariance matrix
    TMatrixDSym *__invCovMatrix;

    int dim=0;
    double determinant; // compute this once and save the cached result
    bool isInverted=false;
    Double_t evaluate() const;
private:
    void invert();
    ClassDef(RooGaussianCorr,1)
};

#endif
