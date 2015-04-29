#ifndef RooParameterizedPdf_hh
#define RooParameterizedPdf_hh

#include "TObject.h"

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsArg.h"

class RooParameterizedPdf : public RooAbsPdf {
public:
  RooParameterizedPdf() {}
  RooParameterizedPdf(const char* name, const char *title,RooRealVar& x, RooAbsReal& scale, RooAbsReal& shift, RooAbsPdf& keysPdf);
  RooParameterizedPdf(const RooParameterizedPdf & other,const char* name=0);

  virtual TObject* clone(const char* newname) const { return new RooParameterizedPdf(*this,newname); }

  inline virtual ~RooParameterizedPdf() {}

protected:
  RooRealProxy __x;
  RooRealProxy __scale;
  RooRealProxy __shift;

  RooRealProxy __pdf;

  Double_t evaluate() const;

private:
  ClassDef(RooParameterizedPdf,1) //no ; on this macro!

    };

#endif
