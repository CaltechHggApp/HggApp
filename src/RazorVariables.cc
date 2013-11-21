#include "RazorVariables.hh"


void RazorVariables::CombineJets(const std::vector<TLorentzVector>& jets,TLorentzVector& hem1, TLorentzVector &hem2) throw(TooManyJets,TooFewJets){
  hem1.SetPtEtaPhiM(0.,0.,0.,0.);
  hem2.SetPtEtaPhiM(0.,0.,0.,0.);
  if(jets.size() > 30) {
    return; //throw new TooManyJets;
  }
  if(jets.size() < 2 ) throw new TooFewJets;

  double M_min=9999999999999.;
  for(int icount=0; icount < pow(2,jets.size()); icount++) { //we will treat icount as a binary array with 0 indicating the jet
    TLorentzVector h1_temp,h2_temp;
    //should be in hem1 and 1 indicating it should be in hem2
    for(int iJet=0;iJet<jets.size();iJet++) {
      if( (icount>>iJet)&1 ) h2_temp+=jets.at(iJet);
      else h1_temp+=jets.at(iJet);
    }
    if(M_min > h1_temp.M2()+h2_temp.M2()) {
      M_min = h1_temp.M2()+h2_temp.M2();
      hem1 = h1_temp;
      hem2 = h2_temp;
    }
  }
}

double RazorVariables::CalcGammaMRstar(const TLorentzVector& hem1, const TLorentzVector& hem2) throw(){
  double A = hem1.P();
  double B = hem2.P();
  double az = hem1.Pz();
  double bz = hem2.Pz();
  TVector3 hem1T, hem2T;
  hem1T.SetXYZ(hem1.Px(),hem1.Py(),0.0);
  hem2T.SetXYZ(hem2.Px(),hem2.Py(),0.0);
  double ATBT = (hem1T+hem2T).Mag2();

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (hem2T.Dot(hem2T)-hem1T.Dot(hem1T))*(hem2T.Dot(hem2T)-hem1T.Dot(hem1T))/(hem1T+hem2T).Mag2());

  double mybeta = (hem2T.Dot(hem2T)-hem1T.Dot(hem1T))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));

  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  //gamma times MRstar                                                                                                                                                                              
  temp *= mygamma;

  return temp;
}


double RazorVariables::CalcMTR(const TLorentzVector& hem1, const TLorentzVector& hem2, const TVector3& met) throw(){

  double temp = met.Mag()*(hem1.Pt()+hem2.Pt()) - met.Dot(hem1.Vect()+hem2.Vect());
  temp /= 2.;

  temp = sqrt(temp);

  return temp;
}
