#include "HggMuonID.hh"

HggMuonID::HggMuonID(){}

bool HggMuonID::passMuMuGammaID(VecbosMu &muon){
  if(!muon.isGlobalMuon || !muon.isTrackerMuon) return false;
  if(muon.nTrackHits <= 10 || muon.nPixelHits==0) return false;
  if(muon.trackImpactPar >=0.2) return false;
  if(muon.trkIso >= 3) return false;
  return true;
}
