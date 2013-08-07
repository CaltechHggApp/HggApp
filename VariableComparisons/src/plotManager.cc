#include "include/plotManager.hh"


plotManager::plotManager(TString inTag):freeze(false),vetos(0),use4Cat(false){
  histNameTag = inTag;
}

void plotManager::addCategory(TString name, TString cut){
  if(freeze) return;
  catNames.push_back(name);
  catCuts.push_back(cut);
}

void plotManager::addVariable(TString name,TString var,int bins, float low, float high,bool correct){
  if(freeze) return;
  varNames.push_back(name);
  variables.push_back(var);
  xBins.push_back(bins);
  xLow.push_back(low);
  xHigh.push_back(high);
  useCorrection.push_back(correct);
}

void plotManager::buildHistograms(){
  if(freeze) return;
  freeze = true;
  for(int i=0;i<catNames.size();i++){
    std::vector<TH1F*> perCat_rp;
    std::vector<TH1F*> perCat_re;
    std::vector<TH1F*> perCat_f;
    for(int j=0;j<variables.size();j++){
      perCat_rp.push_back( new TH1F(varNames.at(j)+"_"+catNames.at(i)+"_"+histNameTag+"_realPho","",xBins.at(j),xLow.at(j),xHigh.at(j)) );
      perCat_re.push_back( new TH1F(varNames.at(j)+"_"+catNames.at(i)+"_"+histNameTag+"_realEle","",xBins.at(j),xLow.at(j),xHigh.at(j)) );
      perCat_f.push_back( new TH1F(varNames.at(j)+"_"+catNames.at(i)+"_"+histNameTag+"_fake","",xBins.at(j),xLow.at(j),xHigh.at(j)) );
      perCat_rp.back()->Sumw2();
      perCat_re.back()->Sumw2();
      perCat_f.back()->Sumw2();
    }
    realPhoHistograms.push_back(perCat_rp);
    realEleHistograms.push_back(perCat_re);
    fakeHistograms.push_back(perCat_f);
  }
}

void plotManager::buildFormulas(TChain *fChain){
  if(!freeze) return;
  etaVal = new TTreeFormula("isEB","abs(etaSC)",(TTree*)fChain);
  isRealPho = new TTreeFormula("isRealPho","realPho==1",(TTree*)fChain);
  isRealEle = new TTreeFormula("isRealEle","realEle==1",(TTree*)fChain);
  isTrigger = new TTreeFormula("isTrigger","Trigger==1",(TTree*)fChain);
  for(int i=0;i<variables.size();i++){
    if(useCorrection.at(i)){
      if(use4Cat) varFormulas.push_back( new TTreeFormula( varNames.at(i)+"_"+histNameTag, corrector4cat.getCorrectString(variables.at(i),  etaVal->EvalInstance()),           (TTree*)fChain ) );
      else        varFormulas.push_back( new TTreeFormula( varNames.at(i)+"_"+histNameTag, corrector.getCorrectString(variables.at(i),      (etaVal->EvalInstance() < 1.48) ), (TTree*)fChain ) );
    }
    else varFormulas.push_back( new TTreeFormula( varNames.at(i)+"_"+histNameTag, variables.at(i), (TTree*)fChain ) );
  }
  for(int i=0;i<catNames.size();i++){
    catFormulas.push_back( new TTreeFormula( catNames.at(i)+"_"+histNameTag, catCuts.at(i), (TTree*)fChain ) );
  }

  if(vetos){
    for(int i=0;i<vetos->size();i++){
      vetoFormulas.push_back(new TTreeFormula(Form("veto%d",i),vetos->at(i),(TTree*)fChain));
    }
  }
}

void plotManager::updateFormulas(){
  std::vector<TTreeFormula*>::iterator it = varFormulas.begin();
  for(; it != varFormulas.end(); it++) (*it)->UpdateFormulaLeaves();
  it = catFormulas.begin();
  for(; it != catFormulas.end(); it++) (*it)->UpdateFormulaLeaves();
  it = vetoFormulas.begin();
  for(; it != vetoFormulas.end(); it++) (*it)->UpdateFormulaLeaves();

  etaVal->UpdateFormulaLeaves();
  isRealPho->UpdateFormulaLeaves();
  isRealEle->UpdateFormulaLeaves();
  isTrigger->UpdateFormulaLeaves();
}

void plotManager::destroyFormulas(){
  varFormulas.clear();
  catFormulas.clear();
  vetoFormulas.clear();
  delete etaVal;
  delete isRealPho;
  delete isRealEle;
  delete isTrigger;
}

void plotManager::processEntry(float weight){
  if(!freeze) return;
  if(! isTrigger->EvalInstance() ) return;
  for(int i=0;i<vetoFormulas.size();i++){
    if( vetoFormulas.at(i)->EvalInstance() ) return;
  }
  for(int i=0;i<catNames.size();i++){    
    if( !catFormulas.at(i)->EvalInstance() ) continue; //not in this category
    for(int j=0;j<variables.size();j++){
      TH1F* hist;
      if(isRealPho->EvalInstance())      hist = realPhoHistograms.at(i).at(j);      
      else if(isRealEle->EvalInstance()) hist = realEleHistograms.at(i).at(j);      
      else                               hist = fakeHistograms.at(i).at(j);      
      hist->Fill(varFormulas.at(j)->EvalInstance(),weight);
    }
  }
}

void plotManager::processChain(TChain *fChain,float weight){
  if(!freeze) buildHistograms();

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
  }

  destroyFormulas();
}

std::array<TH1F*,3> plotManager::getHistogram(TString cat, TString var){
  std::array<TH1F*,3> output = {0,0,0};
  for(int i=0;i<catNames.size();i++){
    if(catNames.at(i)!=cat) continue;
    for(int j=0;j<variables.size();j++){
      if(varNames.at(j) != var) continue;
      output[0] = realPhoHistograms.at(i).at(j);
      output[1] = realEleHistograms.at(i).at(j);
      output[2] = fakeHistograms.at(i).at(j);
      return output;
    }
  }  
  return output;
}

void plotManager::saveAll(TFile *f){
  f->cd();
  for(int i=0;i<catNames.size();i++){
    for(int j=0;j<variables.size();j++){
      realPhoHistograms.at(i).at(j)->Write();
      realEleHistograms.at(i).at(j)->Write();
      fakeHistograms.at(i).at(j)->Write();
    }
  }    
}

