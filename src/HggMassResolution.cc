#include <HggMassResolution.hh>
#include <TLorentzVector.h>
#include "ReadConfig.hh"
#include "TMath.h"

#include <iostream>
using namespace std;
#define debugMassRes 0
HggMassResolution::HggMassResolution(){
  this->clear();

  Categories.push_back("EBlowEtaGoldGAP"); highR9.push_back(true);
  Categories.push_back("EBlowEtaGoldCM"); highR9.push_back(true);
  Categories.push_back("EBlowEtaBad"); highR9.push_back(false);
  Categories.push_back("EBhighEtaGold"); highR9.push_back(true);
  Categories.push_back("EBhighEtaBad"); highR9.push_back(false);
  Categories.push_back("EElowEtaGold"); highR9.push_back(true);
  Categories.push_back("EElowEtaBad"); highR9.push_back(false);
  Categories.push_back("EEhighEtaGold"); highR9.push_back(true);
  Categories.push_back("EEhighEtaBad"); highR9.push_back(false);

  minEta.push_back(0.); maxEta.push_back(1.);
  minEta.push_back(0.); maxEta.push_back(1.);
  minEta.push_back(0.); maxEta.push_back(1.);
  minEta.push_back(1.); maxEta.push_back(1.5);
  minEta.push_back(1.); maxEta.push_back(1.5);
  minEta.push_back(1.5); maxEta.push_back(2.);
  minEta.push_back(1.5); maxEta.push_back(2.);
  minEta.push_back(2.); maxEta.push_back(3.);
  minEta.push_back(2.); maxEta.push_back(3.);

  dzRes.push_back(0.1); dzRes.push_back(TMath::Sqrt(2.)*5.8);

}

void HggMassResolution::clear(){
  
}

double HggMassResolution::getMassResolution(VecbosPho *leadPho,VecbosPho *subleadPho, TVector3 vtx,bool isWrongVtx){
  TLorentzVector p4Pho1 = leadPho->p4FromVtx(vtx,leadPho->finalEnergy);
  TLorentzVector p4Pho2 = subleadPho->p4FromVtx(vtx,subleadPho->finalEnergy);
  double eRes = this->getMassResolutionEonly(leadPho,subleadPho,vtx);
  double angleRes = this->getAngleResolution(leadPho,subleadPho,vtx,isWrongVtx);
  double higgsMass = (p4Pho1+p4Pho2).M();
  
  angleRes*=higgsMass;
  return TMath::Sqrt((eRes*eRes)+(angleRes*angleRes));
}

double HggMassResolution::getMassResolutionEonly(VecbosPho *leadPho,VecbosPho *subleadPho,TVector3 vtx){
  TLorentzVector p4Pho1 = leadPho->p4FromVtx(vtx,leadPho->finalEnergy);
  TLorentzVector p4Pho2 = subleadPho->p4FromVtx(vtx,subleadPho->finalEnergy);
  double resPho1 = this->getResolution(leadPho);
  double resPho2 = this->getResolution(subleadPho);
  double higgsMass = (p4Pho1+p4Pho2).M();
  
  return 0.5*higgsMass*TMath::Sqrt( (resPho1*resPho1)/(leadPho->finalEnergy*leadPho->finalEnergy) + (resPho2*resPho2)/(subleadPho->finalEnergy*subleadPho->finalEnergy) );
}

void HggMassResolution::init(){
  ReadConfig cfg;
  cfg.read(config);

  char valString[400];
  for(int i=0;i<nCategories;i++){
    std::vector<string> thisSmear =  cfg.getTokens(Categories[i],",");
    if(thisSmear.size()!=2){
      cout << "ERROR: INVALID smearing configuration" << endl;
      cout << i << "    " << Categories[i] << "   " << thisSmear[0] << endl;
      throw -100;
      return;
    }
    smear[i] = pair<float,float>(0.01*atof(thisSmear[0].c_str()),0.01*atof(thisSmear[1].c_str()));
  }
}

double HggMassResolution::getAngleResolution(VecbosPho* pho1,VecbosPho* pho2, TVector3 vtx, bool wrongVtx){
  TVector3 pho1Pos= pho1->SC.CaloPos - vtx;
  TVector3 pho2Pos = pho2->SC.CaloPos - vtx;

  double r1 = pho1Pos.Mag();
  double r2 = pho2Pos.Mag();
  double cos = TMath::Cos(pho1Pos.Phi()-pho2Pos.Phi());
  double sech1 = 1./TMath::CosH(pho1Pos.Eta());
  double sech2 = 1./TMath::CosH(pho2Pos.Eta());
  double tanh1 = TMath::TanH(pho1Pos.Eta());
  double tanh2 = TMath::TanH(pho2Pos.Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos);
  double denominator = 1. - tanh1*tanh2 - sech1*sech2*cos;
  if(debugMassRes) cout << ">> >> got Angle Resolution: " << fabs((0.5*dzRes[wrongVtx]/denominator)*(numerator1/r1 + numerator2/r2)) << endl;
  if(debugMassRes) cout << " >> dz: " << dzRes[wrongVtx] << endl;
  return fabs((0.5*dzRes[wrongVtx]/denominator)*(numerator1/r1 + numerator2/r2));
  /*
  //do this exactly like MIT

  Double_t x1= pho1Pos.X();
  Double_t y1= pho1Pos.Y();
  Double_t z1= pho1Pos.Z();

  Double_t x2= pho2Pos.X();
  Double_t y2= pho2Pos.Y();
  Double_t z2= pho2Pos.Z();

  Double_t r1 = sqrt(x1*x1+y1*y1+z1*z1);
  Double_t r2 = sqrt(x2*x2+y2*y2+z2*z2);
  Double_t phi1 = atan2(y1,x1);
  Double_t theta1 = atan2(sqrt(x1*x1+y1*y1),z1);
  Double_t phi2 = atan2(y2,x2);
  Double_t theta2 = atan2(sqrt(x2*x2+y2*y2),z2);

  Double_t sech1 = sin(theta1);
  Double_t tanh1 = cos(theta1);
  Double_t sech2 = sin(theta2);
  Double_t tanh2 = cos(theta2);
  Double_t cos12 = cos(phi1-phi2);

  Double_t rad1 = sech1*(sech1*tanh2-tanh1*sech2*cos12)/(1-tanh1*tanh2-sech1*sech2*cos12);
  Double_t rad2 = sech2*(sech2*tanh1-tanh2*sech1*cos12)/(1-tanh2*tanh1-sech2*sech1*cos12);
  if(debugMassRes) cout << " >> Angle Resolution: " << dzRes[wrongVtx] * 0.5*fabs(rad1/r1 + rad2/r2) << endl;
  return dzRes[wrongVtx] * 0.5*fabs(rad1/r1 + rad2/r2);
  */
}

float HggMassResolution::getResolution(VecbosPho* pho){
  int cat = this->getCategory(pho);
  if(debugMassRes) cout << "category: " << cat << endl;
  pair<float,float> catRes = smear[this->getCategory(pho)];
  if(debugMassRes) cout << "category Res: " << catRes.first << endl;
  return TMath::Sqrt(pho->finalEnergyError*pho->finalEnergyError+
		     catRes.first*catRes.first*pho->finalEnergy*pho->finalEnergy);
}

int HggMassResolution::getCategory(VecbosPho* pho){
  if(debugMassRes) cout << "getCategory" << endl;
  if(this->isSphericalPhoton(pho->SC.BCSeed.iEta,pho->SC.BCSeed.iPhi)){ // these are the "special" photons that have to be treated differently
    if(pho->SC.r9 > r9Cut) return sphericalIndex;
  }
  if(debugMassRes) cout << "not a special photon" << endl;

  for(int i=0;i<nCategories;i++){
    if(i == sphericalIndex) continue;

    if( minEta[i] <= fabs(pho->SC.eta)
	&& fabs(pho->SC.eta) < maxEta[i]
	&& (pho->SC.r9 > r9Cut) == highR9[i]) return i;
  }
  cout << "Photon with no category??? r9: " << pho->SC.r9 << "  eta: " << pho->SC.eta << "  phi: " << pho->SC.phi << endl;
  throw -99;
  return -1;
}

bool HggMassResolution::isSphericalPhoton(int ieta, int iphi)
{                                                             
  if(debugMassRes) cout << "isSphericalPhoton" << endl;
  if ((iphi %20)<=5 || (iphi%20)>=16){                      
    return false;                                         
  }                                                         
  
  int ietaTT=(std::abs(ieta)-1)/5+1;                        
    if                                                        
      (                                                     
       (ietaTT>= 2&&     ietaTT<    5 ) ||               
       (ietaTT>= 7&&     ietaTT<    9 ) ||               
       (ietaTT>= 11&&    ietaTT<    13) ||               
       (ietaTT>= 15&&    ietaTT<    17)                  
       ){                                                
      return true;                                          
    }                                                         
                                                              
                                                              
    return false;                                             
}                                                             
