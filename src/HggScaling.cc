#include "HggScaling.hh"
#include <iostream>
#define debugScaling 0

HggScaling::HggScaling(){
}

HggScaling::HggScaling(std::string cfg){
  this->LoadConfig(cfg);
}

void HggScaling::LoadConfig(std::string cfg){
  ReadConfig cfgReader;
  try{
    cfgReader.read(cfg);
  }catch(const std::runtime_error &e){
    std::cout << "HggScaling" << std::endl;
    throw e;
  }
  if(debugScaling) std::cout << "Loading file " << cfg << std::endl;
  if(debugScaling) cfgReader.printAll();
  std::vector<std::string> allRescales = cfgReader.getAllParameters();

  for(std::vector<std::string>::const_iterator iRescale=allRescales.begin();
      iRescale != allRescales.end(); iRescale++){
    if(debugScaling) std::cout << *iRescale << std::endl;
    std::vector<std::string> s_vals = cfgReader.getTokens(*iRescale,",");
    std::vector<float> vals;
    for(std::vector<std::string>::const_iterator iScale=s_vals.begin();
	iScale != s_vals.end(); iScale++){
      vals.push_back( atof(iScale->c_str()) );
    }
    if(vals.size() != 4) continue; //format needs to be: barrel_lin,barrel_const,endcap_lin,endcap_const
    Corrections[ std::pair<std::string,bool>(*iRescale,0) ] = std::pair<float,float>(vals.at(0),vals.at(1));
    Corrections[ std::pair<std::string,bool>(*iRescale,1) ] = std::pair<float,float>(vals.at(2),vals.at(3));
    std::cout <<"Loading Correction: " << *iRescale
	      << ":  " << vals.at(0) << "  " << vals.at(1)
	      << vals.at(2) << "  " << vals.at(3) << std::endl;
  }
}

void HggScaling::ScalePhoton(VecbosPho &pho){
  CorrectionMap::const_iterator it;
  for(it = Corrections.begin(); it != Corrections.end(); it++){
    ApplyScaling(pho,it->first);
  }
}

void HggScaling::ApplyScaling(VecbosPho &pho, std::pair<std::string,bool> correction){
  if(Corrections.find(correction)==Corrections.end()) return; //don't have a correction of this type
  //get the right correction values
  std::pair<float,float> cor = Corrections[correction];
  std::string corType = correction.first;
    
  if(debugScaling) std::cout << "Applying Correction " << corType << ":   " 
			     << cor.first << "  " << cor.second << std::endl;
  //now we have to be kind of hack-ish:
  if(corType.compare("R9") == 0){
    pho.SC.r9 = cor.first*pho.SC.r9+cor.second;
  }
  if(corType.compare("sieie") == 0){
    pho.SC.sigmaIEtaIEta = cor.first*pho.SC.sigmaIEtaIEta+cor.second;
  }
  if(corType.compare("etaWidth") == 0){
    pho.SC.etaWidth = cor.first*pho.SC.etaWidth+cor.second;
  }
  if(corType.compare("phiWidth") == 0){
    pho.SC.phiWidth = cor.first*pho.SC.phiWidth+cor.second;
  }
  if(corType.compare("s4ratio") == 0){
    pho.SC.s4ratio = cor.first*pho.SC.s4ratio+cor.second;
  }

}
