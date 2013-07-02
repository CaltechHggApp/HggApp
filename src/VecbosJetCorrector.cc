#include "VecbosJetCorrector.hh"

VecbosJetCorrector::VecbosJetCorrector(){}

VecbosJetCorrector::VecbosJetCorrector(std::string cfgStr)
{this->getConfig(cfgStr);}

VecbosJetCorrector::VecbosJetCorrector(ReadConfig& cfg)
{this->getConfig(cfg);}

void VecbosJetCorrector::getConfig(std::string cfgStr){
  ReadConfig cfg(cfgStr);
  this->getConfig(cfg);
}

void VecbosJetCorrector::getConfig(ReadConfig & cfg){
  std::string ResParS = cfg.getParameter("JetCorrection_L2L3Residual");
  std::string L3ParS = cfg.getParameter("JetCorrection_L3Correction");
  std::string L2ParS = cfg.getParameter("JetCorrection_L2Correction");
  std::string L1ParS = cfg.getParameter("JetCorrection_L1Correction");

  vPar.push_back( JetCorrectorParameters(L1ParS) );
  vPar.push_back( JetCorrectorParameters(L2ParS) );
  vPar.push_back( JetCorrectorParameters(L3ParS) );
  vPar.push_back( JetCorrectorParameters(ResParS) );

  jfCorr = new FactorizedJetCorrector(vPar);
}

void VecbosJetCorrector::CorrectJet(VecbosJet& jet,float rho){
  jfCorr->setJetEta(jet.eta);
  jfCorr->setJetPt(TMath::Sqrt(jet.uncorrpx*jet.uncorrpx+jet.uncorrpy*jet.uncorrpy));
  jfCorr->setJetA(jet.area);
  jfCorr->setRho(rho);
  
  double scale = jfCorr->getCorrection();

  jet.energy = jet.uncorrEnergy*scale;

  float px = jet.uncorrpx*scale;
  float py = jet.uncorrpy*scale;

  jet.pt = TMath::Sqrt( px*px + py*py );
}
