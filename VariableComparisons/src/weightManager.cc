#include "include/weightManager.hh"

weightManager::weightManager(){
  crossSection["DiPhotonBox_Pt10to25"] = 424.8;
  crossSection["DiPhotonBox_Pt25to250"] = 15.54;
  crossSection["DiPhotonBox_Pt250"] = 1.18E-3;

  crossSection["DiPhotonJets"] = 75.4;

  crossSection["DiPhotonJets_sherpa"] = 120.354 * 1.1; //k-factor

  crossSection["GJets_Pt20to40"] = 8.19E4*1.84E-3;
  crossSection["GJets_Pt40"] = 8.84E3*5.39E-2;

  crossSection["QCD_Pt30to40"] = 5.2E7*2.35E-4;
  crossSection["QCD_Pt40"] = 2.37E7*2.18E-3;

  crossSection["DYJetsToLL_M-50"] = 3.53E3*2.1; //fake rate correction

  Nevents["DiPhotonBox_Pt10to25"] = 138960;
  Nevents["DiPhotonBox_Pt25to250"] = 493197;
  Nevents["DiPhotonBox_Pt250"] = 464253;

  Nevents["DiPhotonJets"] = 1391248;

  Nevents["DiPhotonJets_sherpa"] = 14403906;
  
  Nevents["GJets_Pt20to40"] = 5605445;
  Nevents["GJets_Pt40"] = 5956139;

  Nevents["QCD_Pt30to40"]=5584126;
  Nevents["QCD_Pt40"]=8996255;  

  Nevents["DYJetsToLL_M-50"] = 30461028;
}

float weightManager::getWeight(TString name,float lumi){
  if(crossSection.find(name)==crossSection.end()) return 0;

  float cs = crossSection[name];
  float N  = Nevents[name];

  float expected = lumi*1000.*cs;
  
  return expected/N;
}

