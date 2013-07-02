// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

// local includes

#include "Vecbos.hh"
#include "ShapeAnalysis.hh"

ShapeAnalysis::ShapeAnalysis(TTree *tree) : Vecbos(tree) {}

ShapeAnalysis::~ShapeAnalysis(){}

vector<TH1D*> ShapeAnalysis::CreateHistos(string dirname){
  vector<TH1D*> histos;
  string name;

 
  return histos;
}

void ShapeAnalysis::FillHistos(vector<TH1D*> histos){

  int ih = -1;
  
 
  
}

void ShapeAnalysis::Loop(string outname) {
  outfilename = outname;
  if(fChain == 0) return;

  vector< vector<TH1D*> > Histos;

  //  for(int i=0; i<41; i++) {
  char name[32];
  //  sprintf(name,"0_0_%i",-20+i);
  sprintf(name,"0_0_%i",0);
  Histos.push_back(CreateHistos(name));
  //}
  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = 100;//fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
    // check the alpgen ID
    
    vector<TLorentzVector> GenCand;
    for(int i = 0; i < nMc; i++){
      if(statusMc[i] == 1){
	if(abs(idMc[i]) != 12 && abs(idMc[i]) != 14 && 
	   abs(idMc[i]) != 16 && abs(idMc[i]) != 13){
	  TLorentzVector J;
	  J.SetPtEtaPhiE(pMc[i]/cosh(etaMc[i]), etaMc[i], phiMc[i], energyMc[i]);
	  GenCand.push_back(J);
	}
      }
    }

    vector<Jet> GenJets = SortJet(SISCone(GenCand, 0.7, 0.0));

    if(nZ0ToMuMu == 0) continue;
    double truemass = 91.1876; // Put the PDG value here
    vector<double> masses;
    for(int iZ =0; iZ<nZ0ToMuMu; iZ++)
      masses.push_back(massZ0ToMuMu[iZ]);
    
    double bestd = 9999999999.;
    int    besti = -999;
    for(int i=0; i<masses.size(); i++) {
      double thisd = fabs(masses[i]-truemass);
      if(thisd < bestd) {
	bestd = thisd;
	besti = i;
      }
    }
    if(besti < 0) continue;
    
    TLorentzVector Z(pxZ0ToMuMu[besti], pyZ0ToMuMu[besti], 
		     pzZ0ToMuMu[besti], energyZ0ToMuMu[besti]);
    TVector3 b = Z.BoostVector();
    b.SetZ(0.0);
    
    //Create Cool Tools object to use functionality
    //when things are calculated (sphericity, thrust, etc.) they 
    //are storing in this object
    CoolTools *CT = new CoolTools();
    
    //Tool boosts jet collection according to velocity vector b (units of c)
    vector<Jet> bj = CT->BoostJets(GenJets, b);

    //Given a vector of type Jet returns the 4Vectors
    //Hence all the tools are written for vectors of 4Vectors
    vector<TLorentzVector> LV = CT->Get4Vectors(bj);
    
    
    
    //Calculates Sphericity - afterwards related values can be accessed
    CT->CalcSphericity(LV);
  
    //Calculates Transverse Sphericity
    CT->CalcTSphericity(LV);

    //Calculates Thrust
    CT->CalcThrust(LV);
 
    //Calculates Tranverse thrust
    CT->CalcTranThrust(LV);
 
    //Calculates Fox-Wolfram variables of all orders
    //Calculates transverse Fox-Wolfram of all orders too
    //Also, spherical harmonic and legendre polynomial functions (see code)
    for(int l = 0; l < 7; l++){
      cout << "Fox-Wolfram order " << l << " = " << CT->FoxWolfram(LV, l) << endl;
      cout << "Tran Fox-Wolfram order " << l << " = " << CT->TranFoxWolfram(LV, l) << endl;
    }

    cout << "Aplanarity " << CT->Aplanarity() << endl;
    cout << "Sphericity " << CT->Sphericity() << endl;
    cout << "TranSphericity " << CT->TranSphericity() << endl;
    cout << "Thrust " << CT->Thrust() << endl;
    cout << "TranThrust " << CT->TranThrust() << endl;
    cout << "ThrustMajor " << CT->ThrustMajor() << endl;
    cout << "ThrustMinor " << CT->ThrustMinor() << endl;
    cout << "TranThrustMinor " << CT->TranThrustMinor() << endl;
    cout << "Oblateness " << CT->Oblateness() << endl;
    cout << "Plus all the axes associated with these things" << endl << endl;

 

  }
    
  TFile *file = new TFile(outfilename.c_str(),"RECREATE");
  //for(int i=0; i<41; i++) {
  //  char name[32];
  sprintf(name,"0_0_%i",0);//-20+i);
  //  WriteHistos(Histos[i], file,name);
  WriteHistos(Histos[0], file,name);
  //} 
}
