#include "MixSpinDatasets.h"
#include <iostream>
MixSpinDatasets::MixSpinDatasets(RooWorkspace *w):ws(w)
{

}

void MixSpinDatasets::scheduleMix(const char* mc1, const char* mc2, 
				  float f1, TString outputName){
  mc1L.push_back(mc1);
  mc2L.push_back(mc2);
  f1L.push_back(f1);
  outputNameL.push_back(outputName);
}

void MixSpinDatasets::mixAll(){
  MakeSpinFits::getLabels("evtcat",&catNames,ws);
  for(int i=0;i<mc1L.size();i++){
    mix(mc1L.at(i).Data(),mc2L.at(i).Data(),f1L.at(i),outputNameL.at(i));
  }
}

void MixSpinDatasets::mix(const char* mc1, const char* mc2, 
			  float f1, TString outputName){
  if(outputName=="") outputName = Form("%s_%0.2f_%s_%0.2f",mc1,f1,mc2,1-f1);
  RooCategory* l = ((RooCategory*)ws->obj("labels"));
  std::cout << outputName << "  " << l << std::endl;
  l->defineType(outputName,l->numBins(""));

  internalMix(mc1,mc2,f1,outputName,"Combined");
  ws->import(*l);
}

void MixSpinDatasets::internalMix(const char* mc1, const char* mc2, 
				  float f1,TString outputName,TString cat){

  RooDataSet *mc1_d = (RooDataSet*)ws->data(Form("%s_%s",mc1,cat.Data()));
  RooDataSet *mc2_d = (RooDataSet*)ws->data(Form("%s_%s",mc2,cat.Data()));

  RooRealVar *evtW = ws->var("evtWeight");

  RooArgSet inputSet(*(mc1_d->get(0)),RooArgSet(*evtW));


  RooDataSet *mixed_TMP = new RooDataSet("TMP","",inputSet); 

  std::cout << "Mixing " << f1 << "*" <<mc1 << "  &&  " << 1-f1 <<"*"<< mc2<<std::endl;

  Long64_t iEntry=-1;
  const RooArgSet *set;

  while( (set=mc1_d->get(++iEntry)) ){ //loop over first DS
    float w = mc1_d->weight();
    float we = mc1_d->weightError();
    evtW->setVal(w*f1);
    evtW->setError(we*f1);
    if(iEntry%5000==0) std::cout << "Adding Entry " << iEntry << "  weight " << evtW->getVal() <<std::endl;
    mixed_TMP->add( RooArgSet(*set,RooArgSet(*evtW)) );
  }

  iEntry=-1;
  while( (set=mc2_d->get(++iEntry)) ){ //loop over first DS
    float w = mc2_d->weight();
    float we = mc2_d->weightError();
    if(iEntry%5000==0) std::cout << "Adding Entry " << iEntry << "  weight " << w*(1-f1) <<std::endl;
    evtW->setVal(w*(1-f1));
    evtW->setError(we*(1-f1));
    mixed_TMP->add( RooArgSet(*set,RooArgSet(*evtW)) );
  }

  RooDataSet *mixed = new RooDataSet(outputName+"_"+cat,"",mixed_TMP,inputSet,0,"evtWeight");

  ws->import( *mixed);

  RooRealVar * tot_EB_1 = ws->var( Form("%s_EB_totalEvents",mc1) );
  RooRealVar * tot_EB_2 = ws->var( Form("%s_EB_totalEvents",mc2) );
  RooRealVar * tot_EE_1 = ws->var( Form("%s_EE_totalEvents",mc1) );
  RooRealVar * tot_EE_2 = ws->var( Form("%s_EE_totalEvents",mc2) );

  RooRealVar * ngen_1   = ws->var( Form("%s_Ngen",mc1) );
  RooRealVar * ngen_2   = ws->var( Form("%s_Ngen",mc2) );
  
  RooRealVar * tot_EB_mix = new RooRealVar( Form("%s_EB_totalEvents",outputName.Data()), "", tot_EB_1->getVal()*f1 + tot_EB_2->getVal()*(1-f1) );
  RooRealVar * tot_EE_mix = new RooRealVar( Form("%s_EE_totalEvents",outputName.Data()), "", tot_EE_1->getVal()*f1 + tot_EE_2->getVal()*(1-f1) );
  RooRealVar * ngen_mix   = new RooRealVar( Form("%s_Ngen",outputName.Data()),           "", ngen_1->getVal()*f1 + ngen_2->getVal()*(1-f1) );

  ws->import(*tot_EB_mix);
  ws->import(*tot_EE_mix);
  ws->import (*ngen_mix);

  delete mixed;
  delete mixed_TMP;
  delete tot_EB_mix;
  delete tot_EE_mix;
  delete ngen_mix;
}
