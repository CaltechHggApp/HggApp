#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooAbsPdf.h"

RooDataHist *subtract(RooRealVar &mass,RooDataHist &dh,const RooArgList &pdf,const RooArgList &frac)
{
    assert(pdf.getSize()==frac.getSize());

    RooArgSet margs(mass);
    RooDataHist *dh_new = dynamic_cast<RooDataHist*>(dh.emptyClone("no_bkg"));
    
    for(int i=0; i<dh.numEntries(); i++ )
    {
        const RooArgSet *args = dh.get(i);
        assert(args->getSize()==1);
        
        const RooAbsArg *x = args->first();
        const RooAbsReal *xx = dynamic_cast<const RooAbsReal *>(x);
        assert(xx!=NULL);
        
        mass = xx->getVal();

        const double bin_size = dh.binVolume();

        double
            weight       = dh.weight(margs,0),
            weight_error = dh.weightError(); //RooAbsData::SumW2);

//        printf("bin=%2d   mass=%g   weight = %g +/- %g\n",i,xx->getVal(),weight,weight_error);

        double v=0;
        for( int j=0; j<pdf.getSize(); j++ )
        {
            const RooAbsPdf  *f = dynamic_cast<const RooAbsPdf *>(&pdf [j]);
            const RooRealVar *c = dynamic_cast<const RooRealVar*>(&frac[j]);
            
//            printf("expected signal (%s): pdf=%g  coef=%g\n",f->GetName(),f->getVal(args),c->getVal());
            
            v += f->getVal(args) * c->getVal() * bin_size;
        }

//        printf("  %g : %g\n",weight,v);
        double new_val = weight-v;
        dh_new->set(*args,new_val,weight_error);
    }
    
    return dh_new;
}

RooDataHist *subtract(RooRealVar &mass, RooDataHist & dh,RooAbsPdf &pdf,float w){
  RooRealVar c("c","",w);
  return subtract(mass,dh,RooArgList(pdf),RooArgList(c));
}
