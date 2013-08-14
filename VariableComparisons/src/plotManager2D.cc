#include "include/plotManager2D.hh"
#include "assert.h"

void plotManager2D::buildHistograms2D(){
  if(freeze2D) return;
  freeze2D = true;
  for(int i=0;i<catNames.size();i++){
    std::vector<std::vector<TH2F*>> perCat_rp;
    std::vector<std::vector<TH2F*>> perCat_re;
    std::vector<std::vector<TH2F*>> perCat_f;
    for(int j=0;j<variables.size();j++){
      std::vector<TH2F*> perCat2D_rp;
      std::vector<TH2F*> perCat2D_re;
      std::vector<TH2F*> perCat2D_f;
      for(int k=0;k<variables.size();k++) {
	if(k>j) {
	  perCat2D_rp.push_back( new TH2F(varNames.at(j)+"_"+varNames.at(k)+"_"+catNames.at(i)+"_"+histNameTag+"_realPho","",xBins.at(j),xLow.at(j),xHigh.at(j),xBins.at(k),xLow.at(k),xHigh.at(k)) );
	  perCat2D_re.push_back( new TH2F(varNames.at(j)+"_"+varNames.at(k)+"_"+catNames.at(i)+"_"+histNameTag+"_realEle","",xBins.at(j),xLow.at(j),xHigh.at(j),xBins.at(k),xLow.at(k),xHigh.at(k)) );
	  perCat2D_f.push_back(  new TH2F(varNames.at(j)+"_"+varNames.at(k)+"_"+catNames.at(i)+"_"+histNameTag+"_fake","",xBins.at(j),xLow.at(j),xHigh.at(j),xBins.at(k),xLow.at(k),xHigh.at(k)) );
	  perCat2D_rp.back()->Sumw2();
	  perCat2D_re.back()->Sumw2();
	  perCat2D_f.back()->Sumw2();
	} else {
	  perCat2D_rp.push_back(0);
	  perCat2D_re.push_back(0);
	  perCat2D_f.push_back(0);
	}
      }
      perCat_rp.push_back(perCat2D_rp);
      perCat_re.push_back(perCat2D_re);
      perCat_f.push_back(perCat2D_f);
    }
    realPho2DHistograms.push_back(perCat_rp);
    realEle2DHistograms.push_back(perCat_re);
    fake2DHistograms.push_back(perCat_f);
  }
}

void plotManager2D::processEntry2D(float weight){
  if(!freeze) return;
  if(! isTrigger->EvalInstance() ) return;
  for(int i=0;i<vetoFormulas.size();i++){
    if( vetoFormulas.at(i)->EvalInstance() ) return;
  }
  for(int i=0;i<catNames.size();i++){    
    if( !catFormulas.at(i)->EvalInstance() ) continue; //not in this category
    for(int j=0;j<variables.size();j++){
      for(int k=j+1;k<variables.size();k++){
	TH2F* hist;
	if(isRealPho->EvalInstance())      hist = realPho2DHistograms.at(i).at(j).at(k);      
	else if(isRealEle->EvalInstance()) hist = realEle2DHistograms.at(i).at(j).at(k);      
	else                               hist = fake2DHistograms.at(i).at(j).at(k);      
	hist->Fill(varFormulas.at(j)->EvalInstance(),varFormulas.at(k)->EvalInstance(),weight);
      }
    }
  }
}

void plotManager2D::processChain(TChain *fChain,float weight){
  if(!freeze) {
    buildHistograms();
    buildHistograms2D();
  }

  buildFormulas(fChain);

  Long64_t iEntry=-1;
  Int_t iTree=-1;
  while( fChain->GetEntry(++iEntry) ){
    if( fChain->GetTreeNumber() != iTree ){
      updateFormulas();
      iTree = fChain->GetTreeNumber();
    }
    if( !(iEntry %10000) ) std::cout << "Processing Entry " << iEntry << std::endl;
    processEntry(weight);
    processEntry2D(weight);
  }

  destroyFormulas();
}


std::array<TH2F*,3> plotManager2D::get2DHistogram(TString cat, TString var1,TString var2){
  std::array<TH2F*,3> output = {0,0,0};
  for(int i=0;i<catNames.size();i++){
    if(catNames.at(i)!=cat) continue;
    for(int j=0;j<variables.size();j++){
      if(varNames.at(j) != var1) continue;
      for(int k=j+1;k<variables.size();k++){
	if(varNames.at(k) != var2) continue;

	int index1=j,index2=k;
	if(realPho2DHistograms.at(i).at(index1).at(index2)==0) std::swap(index1,index2);

	assert(realPho2DHistograms.at(i).at(index1).at(index2)!=0);
	output[0] = realPho2DHistograms.at(i).at(index1).at(index2);
	output[1] = realEle2DHistograms.at(i).at(index1).at(index2);
	output[2] = fake2DHistograms.at(i).at(index1).at(index2);

      }
      return output;
    }
  }  
  return output;
}

void plotManager2D::saveAll2D(TFile *f){
  f->cd();
  for(int i=0;i<catNames.size();i++){
    for(int j=0;j<variables.size();j++){
      for(int k=j+1;k<variables.size();k++) {
	realPho2DHistograms.at(i).at(j).at(k)->Write();
	realEle2DHistograms.at(i).at(j).at(k)->Write();
	fake2DHistograms.at(i).at(j).at(k)->Write();
	TH2F* total = (TH2F*)realPho2DHistograms.at(i).at(j).at(k)->Clone(variables.at(j)+"_"+variables.at(k)+"_"+catNames.at(i)+"_"+histNameTag+"_total");
	total->Add((TH2F*)realEle2DHistograms.at(i).at(j).at(k));
	total->Add((TH2F*)fake2DHistograms.at(i).at(j).at(k));
	total->Write();
      }
    }
  }    
}

