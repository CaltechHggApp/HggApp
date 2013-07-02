#include "HggMCWeight.hh"
//#include "assert.h"
#include "TString.h"
#include <cstdlib>

HggMCWeight::HggMCWeight(){
  //define the allowed weights:
  flags["XSec"]=false;
  flags["BranchingRatio"]=false;
  flags["pileup"]=false;

  functions["XSec"]= &HggMCWeight::getXSec;
  functions["BranchingRatio"]= &HggMCWeight::getBranching;
  functions["pileup"]= &HggMCWeight::getPileup;

}

void HggMCWeight::init(ReadConfig &cfg){

  std::string br = cfg.getParameter("branchingRatios");
  std::string xs = cfg.getParameter("xsecs");
  std::string pu = cfg.getParameter("puWeight");
  branching.read(br);
  xsec.read(xs);

  puWeightFile = new TFile(pu.c_str());
  pileupWeight = (TH1F*)puWeightFile->Get("pileupReWeight");
}

void HggMCWeight::setReweighting(string name, bool val){
  if(flags.find(name) != flags.end())
    flags[name]=val;
}

float HggMCWeight::getWeight(){
  HggMCWeight::flagmap_t::const_iterator flagIt;
  float scale=1;

  //assert(flags.size() == functions.size());

  for(flagIt = flags.begin();
      flagIt != flags.end();
      flagIt++){
    if( !(flagIt->second) ) continue; //skip ignored rescalings
    scale *= (this->*(functions[flagIt->first]))();
  }
  return scale;
}

float HggMCWeight::getXSec(){
  return atof(xsec.getParameter( Form("m%.0f",higgsMass) ).c_str() );
}

float HggMCWeight::getBranching(){
  return atof(branching.getParameter( Form("m%.0f",higgsMass) ).c_str() );
}

float HggMCWeight::getPileup(){
  return pileupWeight->GetBinContent( pileupWeight->FindFixBin(thisPileup) ); 
}
