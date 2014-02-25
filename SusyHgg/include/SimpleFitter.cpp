#include "SimpleFitter.hpp"

SimpleFitter::buildHistograms(int catIndex) {
  
}

TH2F* SimpleFitter::build_histogram(TString name, int catIndex, region reg, float weight) {
  RooRealVar* MR = ws->var("MR");
  RooRealVar* R  = ws->var("Rsq");

  MR->setRange(0,2000.);
  MR->setBins(50);

  R->setRange(0,1.0);
  R->setBins(25);

  TString fullname = Form("%s_%s",name.Data(),cats.at(catIndex).Data());

  TString selectionString = "";

  fitInfo *info = & (per_cat_fit_ranges.at(catIndex));

  switch(reg) {
  case kSignal:
    fullname+="_Signal";
    selectionString = Form("mgg > %f && mgg < %f",info->signal_min,info->signal_max);
    break;
  case kSideband:
    fullname+="_Sideband";
    selectionString = Form("(mgg > %f && mgg < %f) || (mgg > %f && mgg < %f)",
			   info->sideband_low_min,  info->sideband_low_max,
			   info->sideband_high_min, info->sideband_high_max);
    break;
  case kAll:
    fullname+="_All";
    selectionString = "1";
    break;
  default:
    throw new std::runtime_error( Form("invalid region %d",reg) );
  }

  RooDataHist roohist(Form("%s_hist",fullname.Data()),"",RooArgSet(*MR,*R),*((RooDataSet*)ws->data(Form("%s_%s",name.Data(),cats.at(catIndex).Data()))->reduce(selectionString)),
		   weight);
  RooHistPdf roopdf(Form("%s_pdf",fullname.Data()),"",RooArgSet(*MR,*R),roohist);

  std::cout << roohist.GetName() << std::endl;

  ws->import(roohist);
  ws->import(roopdf);


  return hist;

}
