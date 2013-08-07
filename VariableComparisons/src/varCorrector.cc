#include "include/varCorrector.hh"

varCorrector::varCorrector(){
  correctionsEB["r9"] = "1.0045*r9+0.001";
  correctionsEE["r9"] = "1.0086*r9-0.0007";

   correctionsEB["sieie"] = "0.891832*sieie+0.0009133";
   correctionsEE["sieie"] = "0.9947*sieie  +0.00003";

  correctionsEB["sipip"] = "0.993*sipip-0.0";

  correctionsEB["etaWidth"] = "1.04302 *etaWidth-0.000618";
  correctionsEE["etaWidth"] = "0.903254*etaWidth+0.001346";

  correctionsEB["phiWidth"] = "1.00002*phiWidth-0.000371";
  correctionsEE["phiWidth"] = "0.99992*phiWidth-4.8e-7";

  correctionsEB["e3x3/energyBC"] = "("+correctionsEB["r9"]+")*rawE/energyBC";
  correctionsEE["e3x3/energyBC"] = "("+correctionsEE["r9"]+")*rawE/energyBC";

  correctionsEB["e5x5/energyBC"] = "min(1.0,e5x5/energyBC*1.0022)";

  correctionsEB["e5x5/rawE"] = "("+correctionsEB["e5x5/energyBC"]+")*energyBC/rawE";
  correctionsEE["e5x5/rawE"] = "min(1.0,1.022*e5x5/rawE)";

  correctionsEB["eMax/energyBC"] = "1.012*eMax/energyBC";
  correctionsEE["eMax/energyBC"] = "1.005*eMax/energyBC";
  
  correctionsEB["e2nd/energyBC"] = "1.00*e2nd/energyBC";
  correctionsEE["e2nd/energyBC"] = "1.02*e2nd/energyBC";
  
  correctionsEB["eTop/energyBC"] = "0.94*eTop/energyBC";
  correctionsEE["eTop/energyBC"] = "0.96*eTop/energyBC";
  
  correctionsEB["eBottom/energyBC"] = "0.94*eBottom/energyBC";
  correctionsEE["eBottom/energyBC"] = "0.96*eBottom/energyBC";
  
  correctionsEB["eLeft/energyBC"] = "0.94*eLeft/energyBC";
  correctionsEE["eLeft/energyBC"] = "0.96*eLeft/energyBC";
  
  correctionsEB["eRight/energyBC"] = "0.94*eRight/energyBC";
  correctionsEE["eRight/energyBC"] = "0.96*eRight/energyBC";
  
  correctionsEB["e2x5Max/energyBC"] = "1.006 *e2x5Max/energyBC";
  correctionsEE["e2x5Max/energyBC"] = "1.0075*e2x5Max/energyBC";

  correctionsEB["e2x5Top/energyBC"] = "1.09*e2x5Top/energyBC";
  correctionsEE["e2x5Top/energyBC"] = "1.13*e2x5Top/energyBC";

  correctionsEB["e2x5Bottom/energyBC"] = "1.09*e2x5Bottom/energyBC";
  correctionsEE["e2x5Bottom/energyBC"] = "1.13*e2x5Bottom/energyBC";

  correctionsEB["e2x5Left/energyBC"] = "1.09*e2x5Left/energyBC";
  correctionsEE["e2x5Left/energyBC"] = "1.13*e2x5Left/energyBC";

  correctionsEB["e2x5Right/energyBC"] = "1.09*e2x5Right/energyBC";
  correctionsEE["e2x5Right/energyBC"] = "1.13*e2x5Right/energyBC";

  correctionsEE["se"] = "(etaSC<1.8)*(0.1845*se*se+se-0.00)+(etaSC>=1.8)*(0.168*se*se+se)";
  correctionsEB["se"] = "(etaSC<0.8)*(0.11*se*se+1.025*se)+(etaSC>=0.8)*(0.044*se*se+1.025*se)";
}

TString varCorrector::getCorrectString(TString var, bool isEB){
  std::map<TString,TString> *themap;
  if( isEB ) themap = &correctionsEB;
  else themap = &correctionsEE;

  if( themap->find(var) == themap->end() ) return var; // no correction
  return (*themap)[var];
}




varCorrector4cat::varCorrector4cat(){
  correctionsEBlow["sieie"] = "0.991833*sieie+1.55e-06";
  correctionsEBhigh["sieie"] = "0.999988*sieie-5.777e-5";
  correctionsEElow["sieie"] = "sieie";
  correctionsEElow["sieie"] = "0.9912*sieie+1.4477e-8";

  correctionsEBlow["sieip"] = "0.975*sieip";
  correctionsEBhigh["sieip"] = "0.95*sieip";
  correctionsEElow["sieip"] = "1.025*sieip";
  correctionsEElow["sieip"] = "1.05*sieip";

  correctionsEBlow["sipip"] = "sipip";
  correctionsEBhigh["sipip"] = "1.00162*sipip";
  correctionsEElow["sipip"] = "0.975*sipip";
  correctionsEElow["sipip"] = "1.025*sipip";

  correctionsEBlow["etaWidth"] = "0.975*etaWidth";
  correctionsEBhigh["etaWidth"] = "0.975*etaWidth";
  correctionsEElow["etaWidth"] = "0.947114*etaWidth+0.00020317";
  correctionsEElow["etaWidth"] = "etaWidth+0.0015";

  correctionsEBlow["phiWidth"] = "1.0043*phiWidth-2.79585e-5";
  correctionsEBhigh["phiWidth"] = "0.999763*phiWidth-3.76816e-05";
  correctionsEElow["phiWidth"] = "0.975*phiWidth-0.0002414";
  correctionsEElow["phiWidth"] = "1.00904*phiWidth+0.00182084";

  correctionsEBlow["HE"] = "0.975*HE";
  correctionsEBhigh["HE"] = "0.975*HE+7.6055e-5";
  correctionsEElow["HE"] = "0.85*HE";
  correctionsEElow["HE"] = "1.025*HE-0.0015";

  correctionsEBlow["r9"] = "1*r9+0.0045";
  correctionsEBhigh["r9"] = "0.993464*r9+0.0075";
  correctionsEElow["r9"] = "1.00381*r9-0.0165";
  correctionsEElow["r9"] = "1*r90.006";
}

TString varCorrector4cat::getCorrectString(TString var, float eta){
  std::map<TString,TString> *themap;
  if( fabs(eta) < 1 ) themap = &correctionsEBlow;
  else if( fabs(eta) < 1.44 ) themap = &correctionsEBhigh;
  else if( fabs(eta) < 2. ) themap = &correctionsEElow;
  else themap = &correctionsEEhigh;

  if( themap->find(var) == themap->end() ) return var; // no correction
  return (*themap)[var];
}
