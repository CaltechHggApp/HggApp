#include "TTreeFormula.h"
#include "TString.h"
#include "TChain.h"

typedef std::map<TString,TTreeFormula*> CatMap;

void updateFormulas(CatMap* m){
  CatMap::iterator it;
  for(it=m->begin(); it!=m->end(); it++) (*it).second->UpdateFormulaLeaves();
}

int getCat(CatMap* cats, TString type, const int nCat){
  if((*cats).find(Form("%s_R",type.Data())) != (*cats).end())  if( (*cats)[ Form("%s_R",type.Data()) ]->EvalInstance() ){
    return -1; // rejections
  }
  if((*cats).find(Form("%s_PreSelEta",type.Data())) != (*cats).end()) if( (*cats)[ Form("%s_PreSelEta",type.Data()) ]->EvalInstance() ){
    return -1; // rejection photons between EB and EE
  }
  if((*cats).find(Form("%s_PreSelEt",type.Data())) != (*cats).end()) if( (*cats)[ Form("%s_PreSelEt",type.Data()) ]->EvalInstance() ){
    return -1; // final pT requirements on the photons
  }
  

  for(int iCat=0;iCat<nCat;iCat++){
    if( (*cats)[ Form("%s_%d",type.Data(),iCat) ]->EvalInstance() ) return iCat;
  }
  return -1;
}

CatMap getCategoryCuts(TChain* fChain){
  std::map<TString,TTreeFormula*> categories;


  //categories["multiplicity"] = new TTreeFormula("Mult","Photon>=2",fChain);
  //categories["MVA_PreSelEta"] = new TTreeFormula("MVA_PreSelEta","(abs(Photon[0].etaSC) > 1.4442 && abs(Photon[0].etaSC) < 1.566) || (abs(Photon[1].etaSC) > 1.4442 && abs(Photon[1].etaSC) < 1.566) || abs(Photon[0].etaSC) > 2.5 || abs(Photon[1].etaSC) > 2.5",fChain);
  //categories["MVA_PreSelEt"] = new TTreeFormula("MVA_PreSelEt","Photon[0].pt < 25. || Photon[1].pt < 25. || (Photon[0].pt < 33.3 && Photon[1].pt < 33.3) || Photon[0].pt/mPair < 30./120. || Photon[1].pt/mPair < 40./120.",fChain);
  //categories["MVA_PreSelEt"] = new TTreeFormula("MVA_PreSelEt","(Photon[0].pt/mPair < 30/120) || (Photon[1].pt/mPair < 30/120) || (Photon[0].pt/mPair < 40./120. && Photon[1].pt/mPair < 40./120.)",fChain);
  categories["MVA_R"] = new TTreeFormula("MVA_R","mPair<0 || diPhotonMVA<-0.05 || nPhoton!=2",fChain); // rejection criterion
  categories["MVA_0"] = new TTreeFormula("MVA_0"," (Mjj >= 500) &&  (diPhotonMVA>=-0.05) && (ptJet1 > 30. && ptJet2 > 30.)",fChain);  
  categories["MVA_1"] = new TTreeFormula("MVA_1"," (Mjj >= 250) &&  (diPhotonMVA>=-0.05)",fChain);  
  categories["MVA_2"] = new TTreeFormula("MVA_2"," (diPhotonMVA>=0.88)",fChain);  
  categories["MVA_3"] = new TTreeFormula("MVA_3"," (diPhotonMVA>=0.71) && (diPhotonMVA<0.88)",fChain);  
  categories["MVA_4"] = new TTreeFormula("MVA_4"," (diPhotonMVA>=0.50) && (diPhotonMVA<0.71)",fChain);  
  categories["MVA_5"] = new TTreeFormula("MVA_5"," (diPhotonMVA>=-0.05) && (diPhotonMVA<0.50)",fChain);  
  
  //categories["PFCiC_PreSelEta"] = new TTreeFormula("PFCiC_PreSelEta","(abs(PhotonPFCiC[0].etaSC) > 1.4442 && abs(PhotonPFCiC[0].etaSC) < 1.566) || (abs(PhotonPFCiC[1].etaSC) > 1.4442 && abs(PhotonPFCiC[1].etaSC) < 1.566) || abs(PhotonPFCiC[0].etaSC) > 2.5 || abs(PhotonPFCiC[1].etaSC) > 2.5",fChain);
  //categories["PFCiC_PreSelEt"] = new TTreeFormula("PFCiC_PreSelEt","PhotonPFCiC[0].pt < 25 || PhotonPFCiC[1].pt < 25 || (PhotonPFCiC[0].pt < 33. && PhotonPFCiC[1].pt < 33.) || PhotonPFCiC[0].pt/mPairPFCiC < 30./120. || PhotonPFCiC[1].pt/mPairPFCiC < 40./120.",fChain);
  //categories["PFCiC_PreSelEt"] = new TTreeFormula("MVA_PreSelEt","(PhotonPFCiC[0].pt/mPairPFCiC < 30/120) || (PhotonPFCiC[1].pt/mPairPFCiC < 30/120) || (PhotonPFCiC[0].pt/mPairPFCiC < 40./120. && PhotonPFCiC[1].pt/mPairPFCiC < 40./120.)",fChain);
  categories["PFCiC_R"] = new TTreeFormula("PFCiC_R","mPairPFCiC<0 || nPhotonPFCiC!=2",fChain); // rejection criterion
  categories["PFCiC_0"] = new TTreeFormula("PFCiC_0","(MjjPFCiC >= 500) && (ptJet1 > 30.&& ptJet2 > 30.)",fChain);
  categories["PFCiC_1"] = new TTreeFormula("PFCiC_1","(MjjPFCiC >= 250)",fChain);
  categories["PFCiC_2"] = new TTreeFormula("PFCiC_2","(PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && (abs(PhotonPFCiC[0].etaSC) < 1.48 && abs(PhotonPFCiC[1].etaSC) < 1.48)",fChain);
  categories["PFCiC_3"] = new TTreeFormula("PFCiC_3","!(PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && (abs(PhotonPFCiC[0].etaSC) < 1.48 && abs(PhotonPFCiC[1].etaSC) < 1.48)",fChain);
  categories["PFCiC_4"] = new TTreeFormula("PFCiC_4","(PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && !(abs(PhotonPFCiC[0].etaSC) < 1.48 && abs(PhotonPFCiC[1].etaSC) < 1.48)",fChain);
  categories["PFCiC_5"] = new TTreeFormula("PFCiC_5","!(PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && !(abs(PhotonPFCiC[0].etaSC) < 1.48 && abs(PhotonPFCiC[1].etaSC) < 1.48)",fChain);

  //categories["CiC_PreSelEta"] = new TTreeFormula("CiC_PreSelEta","(abs(PhotonCiC[0].etaSC) > 1.4442 && abs(PhotonCiC[0].etaSC) < 1.566) || (abs(PhotonCiC[1].etaSC) > 1.4442 && abs(PhotonCiC[1].etaSC) < 1.566) || abs(PhotonCiC[0].etaSC) > 2.5 || abs(PhotonCiC[1].etaSC) > 2.5",fChain);
  //categories["CiC_PreSelEt"] = new TTreeFormula("CiC_PreSelEt","PhotonCiC[0].pt < 25 || PhotonCiC[1].pt < 25 || (PhotonCiC[0].pt < 33 && PhotonCiC[1].pt < 33)",fChain);
  categories["CiC_R"] = new TTreeFormula("CiC_R","mPairCiC==-1 || nPhotonCiC==0",fChain); // rejection criterion
  categories["CiC_0"] = new TTreeFormula("CiC_0","(MjjCiC >= 500) && (PhotonCiC[0].pt>30. && PhotonCiC[1].pt >30.) && (ptJet1 > 30. && ptJet2 > 30.)",fChain);
  categories["CiC_1"] = new TTreeFormula("CiC_1","(MjjCiC >= 250) && (PhotonCiC[0].pt>30. || PhotonCiC[1].pt >30.)",fChain);
  categories["CiC_2"] = new TTreeFormula("CiC_2","(PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && (abs(PhotonCiC[0].etaSC) < 1.48 && abs(PhotonCiC[1].etaSC) < 1.48)",fChain);
  categories["CiC_3"] = new TTreeFormula("CiC_3","!(PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && (abs(PhotonCiC[0].etaSC) < 1.48 && abs(PhotonCiC[1].etaSC) < 1.48)",fChain);
  categories["CiC_4"] = new TTreeFormula("CiC_4","(PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && !(abs(PhotonCiC[0].etaSC) < 1.48 && abs(PhotonCiC[1].etaSC) < 1.48)",fChain);
  categories["CiC_5"] = new TTreeFormula("CiC_5","!(PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && !(abs(PhotonCiC[0].etaSC) < 1.48 && abs(PhotonCiC[1].etaSC) < 1.48)",fChain);  

  CatMap::iterator it;
  for(it=categories.begin(); it!=categories.end(); it++) (*it).second->SetQuickLoad(true);  


  return categories;
}
