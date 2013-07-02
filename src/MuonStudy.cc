// std includes
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>

using namespace std;

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <THnSparse.h>

// VecbosApp includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "MuonStudy.hh"
#include "Jet.hh"
#include "CaloTower.hh"

MuonStudy::MuonStudy(TTree *tree, string outname) : Vecbos(tree) {
  outfilename = outname;
}

MuonStudy::~MuonStudy(){}

void MuonStudy::WriteHistosND(vector<THnSparseD*> histos, TFile* file, string dirname){
  file->cd();
  if(!file->GetDirectory(dirname.c_str()))
    file->mkdir(dirname.c_str());
  file->cd(dirname.c_str());
  for(int i=0; i< int(histos.size()); i++) 
    histos[i]->Write();
  file->cd();
}

vector<THnSparseD*> MuonStudy::CreateHistosND(string dirname){
  vector<THnSparseD*> histos;

  int theBins[5];  
  double theMin[5];
  double theMax[5];

  theBins[0] = 50; // sumPt05
  theBins[1] = 50; // emEt05
  theBins[2] = 50; // hadEt05
  theBins[3] = 50; // |dxySignificance|
  theBins[4] = 50; // |dszPV|
  // 
  theMin[0] =  0.;
  theMin[1] =  0.; 
  theMin[2] =  0.; 
  theMin[3] =  0.; 
  theMin[4] =  0.; 
  //
  theMax[0] = 160.;
  theMax[1] = 160.;
  theMax[2] = 180.;
  theMax[3] = 20.;
  theMax[4] = 10.;

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    string name = dirname+"H_NVarDim_"+string(numjets);
    histos.push_back(new THnSparseD(name.c_str(),name.c_str(), 5, theBins, theMin, theMax));
  }

  return histos;
}

vector<TH2D*> MuonStudy::CreateHistos2D(string dirname){
  vector<TH2D*> histos;
  string name;
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt03_vs_MET_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 200., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt03_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt03_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt03_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05_vs_MET_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 200., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt05_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
    for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05-sumPt03_vs_MET_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 200., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05-emEt03_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05-hadEt03_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt05-hoEt03_vs_MET__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }

  ////////////////////////////////////////////

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt03_vs_mT_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 200., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt03_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt03_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt03_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05_vs_mT_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 200., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt05_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
    for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05-sumPt03_vs_mT_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 200., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05-emEt03_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05-hadEt03_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt05-hoEt03_vs_mT__"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 200., 100, 0., 80.));
  }

  ////////////////////////////////////////////////

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt03_vs_sumPt03_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt03_vs_sumPt03_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt03_vs_emEt03_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05_vs_sumPt05_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05_vs_sumPt05_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05_vs_emEt05_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05-emEt03_vs_sumPt05-sumPt03_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05-hadEt03_vs_sumPt05-sumPt03_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05-hadEt03_vs_emEt05-emEt03_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, 0., 80., 100, 0., 80.));
  }

  // Cross correlation

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05_vs_dxySignificance_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, -20., 20., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05_vs_dxySignificance_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, -20., 20., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05_vs_dxySignificance_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, -20., 20., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05_vs_dszPV_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, -1., 1., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05_vs_dszPV_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, -1., 1., 100, 0., 80.));
  }
  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05_vs_dszPV_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, -1., 1., 100, 0., 80.));
  }
  
  ///////////////////////////////////////////////////////////////////////////////

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dszPV_vs_dxySignificance_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(),  100, -20., 20., 100, -1., 1.));
  }
  
  return histos;

}

vector<TH1D*> MuonStudy::CreateHistos1D(string dirname){
  vector<TH1D*> histos;
  string name;

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dxyPV_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -1.0, 1.0));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dxyPVSignificance_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -20., 20.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dszPV_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 300, -1., 1.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dszPVSignificance_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -20., 20.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dxy_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -1.0, 1.0));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dxySignificance_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -20., 20.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dsz_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 300, -1., 1.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dszSignificance_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -20., 20.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_DeltaX_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -0.1, 0.1));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_DeltaY_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -0.1, 0.1));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_DeltaZ_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 300, -15., 15.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_DeltaZSignificance_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 300, -15., 15.));
  }

  // PT

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_pT_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 300.));
  }

  // ISOLATION

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 20.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_nTrk03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 15, 0., 15.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt05_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 20.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_nTrk05_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 15, 0., 15.));
  }



  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_sumPt05-03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_emEt05-03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hadEt05-03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 80.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_hoEt05-03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 20.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_nTrk05-03_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 15, 0., 15.));
  }



  // closest jet

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_dR_Mu_jet_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 50, 0., 5.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_JetIndex_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 11, -0.5, 10.5));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_jetpT_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 300.));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_DeltaPhi_Mu_Jet_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 2.*asin(1.)));
  }

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_nMuon_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 5, -0.5, 4.5));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = dirname+"_ptMuon_"+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0.0, 100.0));
  }
    
  name = dirname+"_MET";
  histos.push_back(new TH1D(name.c_str(),name.c_str(),100, 0., 200.));

  name = dirname+"_METphi";
  histos.push_back(new TH1D(name.c_str(),name.c_str(),100, 0., asin(1)*2.));
  
  name = dirname+"_nJets";
  histos.push_back(new TH1D(name.c_str(),name.c_str(),6, 0., 6.));

  return histos;
}

void MuonStudy::FillHistos(vector<TH1D*> histos, vector<TH2D*> histos2d, vector<THnSparseD*> histosNd, int iMu, vector<Jet> jets ){

  njets = jets.size();

  for(int ij=1; ij<6; ij++) {
    if(njets>= ij) {
      // Find the Track of the Muons
      int iTk = FindTrackMu(iMu);
      if(iTk != -99) {// The track was found
	histos[135+(ij-1)]->Fill(dRMin,weight);

	if(nPV>0) { // A PV was found

	  // 1D Histograms wrt highest SumpT PV
	  histos[0+(ij-1)]->Fill(trackDxyPVTrack[iTk],weight); // signed impact parameter
	  histos[5+(ij-1)]->Fill(trackDxyPVTrack[iTk]/trackDxyErrorTrack[iTk],weight); // signed significance of the impact parameter
	  histos[10+(ij-1)]->Fill(trackDszPVTrack[iTk],weight); // signed SZ distance 
	  histos[15+(ij-1)]->Fill(trackDszPVTrack[iTk]/trackDszErrorTrack[iTk],weight); // signed significance of  signed SZ distance
	  
	  // 1D Histograms wrt Origin
	  histos[20+(ij-1)]->Fill(trackDxyTrack[iTk],weight); // signed impact parameter
	  histos[25+(ij-1)]->Fill(trackDxyTrack[iTk]/trackDxyErrorTrack[iTk],weight); // signed significance of the impact parameter
	  histos[30+(ij-1)]->Fill(trackDszTrack[iTk],weight); // signed SZ distance
	  histos[35+(ij-1)]->Fill(trackDszTrack[iTk]/trackDszErrorTrack[iTk],weight); // signed significance of  signed SZ distance
	  
	  // 1D Histogram from (0,0,0)
	  histos[40+(ij-1)]->Fill(vertexXTrack[iTk]-PVxPV[0],weight);
	  histos[45+(ij-1)]->Fill(vertexYTrack[iTk]-PVyPV[0],weight);
	  histos[50+(ij-1)]->Fill(vertexZTrack[iTk]-PVzPV[0],weight);
	  histos[55+(ij-1)]->Fill((vertexXTrack[iTk]-PVxPV[0])/trackDzErrorTrack[iTk],weight);
	  
	  // pT plots
	  histos[60+(ij-1)]->Fill(pT(pxTrack[iTk],pyTrack[iTk]),weight); // pT
	  histos[165+(ij-1)]->Fill(pT(pxMuon[iMu],pyMuon[iMu]),weight);

	  histos2d[165+(ij-1)]->Fill(trackDxyPVTrack[iTk]/trackDxyErrorTrack[iTk],emEt05Muon[iMu],weight);
	  histos2d[170+(ij-1)]->Fill(trackDxyPVTrack[iTk]/trackDxyErrorTrack[iTk],hadEt05Muon[iMu],weight);
	  histos2d[175+(ij-1)]->Fill(trackDxyPVTrack[iTk]/trackDxyErrorTrack[iTk],sumPt05Muon[iMu],weight);
	  histos2d[180+(ij-1)]->Fill(trackDszTrack[iTk],emEt05Muon[iMu],weight);
	  histos2d[185+(ij-1)]->Fill(trackDszTrack[iTk],hadEt05Muon[iMu],weight);
	  histos2d[190+(ij-1)]->Fill(trackDszTrack[iTk],sumPt05Muon[iMu],weight);
	  histos2d[195+(ij-1)]->Fill(trackDxyPVTrack[iTk]/trackDxyErrorTrack[iTk],trackDszTrack[iTk],weight);

	  double toFill[5];
	  toFill[0] = sumPt05Muon[iMu];
	  toFill[1] = emEt05Muon[iMu];
	  toFill[2] = hadEt05Muon[iMu];
	  toFill[3] = trackDxyTrack[iTk]/trackDxyErrorTrack[iTk];
	  toFill[4] = trackDszPVTrack[iTk];
	  histosNd[ij-1]->Fill(toFill,weight);

	}
      }

      histos[65+(ij-1)]->Fill(sumPt03Muon[iMu],weight);
      histos[70+(ij-1)]->Fill(emEt03Muon[iMu],weight);
      histos[75+(ij-1)]->Fill(hadEt03Muon[iMu],weight);
      histos[80+(ij-1)]->Fill(hoEt03Muon[iMu],weight);
      histos[85+(ij-1)]->Fill(nTrk03Muon[iMu],weight);
      
      histos[90+(ij-1)]->Fill(sumPt05Muon[iMu],weight);
      histos[95+(ij-1)]->Fill(emEt05Muon[iMu],weight);
      histos[100+(ij-1)]->Fill(hadEt05Muon[iMu],weight);
      histos[105+(ij-1)]->Fill(hoEt05Muon[iMu],weight);
      histos[110+(ij-1)]->Fill(nTrk05Muon[iMu],weight);
      
      histos[115+(ij-1)]->Fill(sumPt05Muon[iMu]-sumPt03Muon[iMu],weight);
      histos[120+(ij-1)]->Fill(emEt05Muon[iMu]-emEt03Muon[iMu],weight);
      histos[125+(ij-1)]->Fill(hadEt05Muon[iMu]-hadEt03Muon[iMu],weight);
      histos[130+(ij-1)]->Fill(hoEt05Muon[iMu]-hoEt03Muon[iMu],weight);
      histos[135+(ij-1)]->Fill(nTrk05Muon[iMu]-nTrk03Muon[iMu],weight);
      
      // MET
      histos2d[(ij-1)]->Fill(etMet[0],sumPt03Muon[iMu],weight);
      histos2d[5+(ij-1)]->Fill(etMet[0],emEt03Muon[iMu],weight);
      histos2d[10+(ij-1)]->Fill(etMet[0],hadEt03Muon[iMu],weight);
      histos2d[15+(ij-1)]->Fill(etMet[0],hoEt03Muon[iMu],weight);
      
      histos2d[20+(ij-1)]->Fill(etMet[0],sumPt05Muon[iMu],weight);
      histos2d[25+(ij-1)]->Fill(etMet[0],emEt05Muon[iMu],weight);
      histos2d[30+(ij-1)]->Fill(etMet[0],hadEt05Muon[iMu],weight);
      histos2d[35+(ij-1)]->Fill(etMet[0],hoEt05Muon[iMu],weight);
      
      histos2d[40+(ij-1)]->Fill(etMet[0],sumPt05Muon[iMu]-sumPt03Muon[iMu],weight);
      histos2d[45+(ij-1)]->Fill(etMet[0],emEt05Muon[iMu]-emEt03Muon[iMu],weight);
      histos2d[50+(ij-1)]->Fill(etMet[0],hadEt05Muon[iMu]-hadEt03Muon[iMu],weight);
      histos2d[55+(ij-1)]->Fill(etMet[0],hoEt05Muon[iMu]-hoEt03Muon[iMu],weight);
      
      // mT
      histos2d[60+(ij-1)]->Fill(mT(iMu),sumPt03Muon[iMu],weight);
      histos2d[65+(ij-1)]->Fill(mT(iMu),emEt03Muon[iMu],weight);
      histos2d[70+(ij-1)]->Fill(mT(iMu),hadEt03Muon[iMu],weight);
      histos2d[75+(ij-1)]->Fill(mT(iMu),hoEt03Muon[iMu],weight);
      
      histos2d[80+(ij-1)]->Fill(mT(iMu),sumPt05Muon[iMu],weight);
      histos2d[85+(ij-1)]->Fill(mT(iMu),emEt05Muon[iMu],weight);
      histos2d[90+(ij-1)]->Fill(mT(iMu),hadEt05Muon[iMu],weight);
      histos2d[95+(ij-1)]->Fill(mT(iMu),hoEt05Muon[iMu],weight);
      
      histos2d[100+(ij-1)]->Fill(mT(iMu),sumPt05Muon[iMu]-sumPt03Muon[iMu],weight);
      histos2d[105+(ij-1)]->Fill(mT(iMu),emEt05Muon[iMu]-emEt03Muon[iMu],weight);
      histos2d[110+(ij-1)]->Fill(mT(iMu),hadEt05Muon[iMu]-hadEt03Muon[iMu],weight);
      histos2d[115+(ij-1)]->Fill(mT(iMu),hoEt05Muon[iMu]-hoEt03Muon[iMu],weight);
      
      // correlation among variables
      
      histos2d[120+(ij-1)]->Fill(sumPt03Muon[iMu],emEt03Muon[iMu],weight);
      histos2d[125+(ij-1)]->Fill(sumPt03Muon[iMu],hadEt03Muon[iMu],weight);
      histos2d[130+(ij-1)]->Fill(emEt03Muon[iMu],hadEt03Muon[iMu],weight);
      
      histos2d[135+(ij-1)]->Fill(sumPt05Muon[iMu],emEt05Muon[iMu],weight);
      histos2d[140+(ij-1)]->Fill(sumPt05Muon[iMu],hadEt05Muon[iMu],weight);
      histos2d[145+(ij-1)]->Fill(emEt05Muon[iMu],hadEt05Muon[iMu],weight);
      
      histos2d[150+(ij-1)]->Fill(sumPt05Muon[iMu]-sumPt03Muon[iMu],emEt05Muon[iMu]-emEt03Muon[iMu],weight);
      histos2d[155+(ij-1)]->Fill(sumPt05Muon[iMu]-sumPt03Muon[iMu],hadEt05Muon[iMu]-hadEt03Muon[iMu],weight);
      histos2d[160+(ij-1)]->Fill(emEt05Muon[iMu]-emEt03Muon[iMu],hadEt05Muon[iMu]-hadEt03Muon[iMu],weight);

      // find the closest jet passing the cuts
      double dRmin_mu_jet = 99999999999.;
      int iJetMin = -99;
      for(int k=0; k< njets; k++ ) {
	double dR_mu_jet = DeltaR(double(etaMuon[iMu]),double(phiMuon[iMu]),jets[k].eta(),jets[k].phi()); 
	if(dR_mu_jet<dRmin_mu_jet) {
	  dRmin_mu_jet = dR_mu_jet;
	  iJetMin = k;
	}
      }
      
      if(iJetMin != -99) {
	histos[140+(ij-1)]->Fill(dRmin_mu_jet,weight);
	histos[145+(ij-1)]->Fill(iJetMin,weight);
	histos[150+(ij-1)]->Fill(jets[iJetMin].pt(),weight);
	double dphi = DeltaPhi(double(phiMuon[iMu]),jets[iJetMin].phi());
	histos[155+(ij-1)]->Fill(fabs(dphi),weight);
      }
      
      histos[160+(ij-1)]->Fill(nMuon,weight);
    }
  }

  histos[int(histos.size()-3)]->Fill(etMet[0],weight);
  histos[int(histos.size()-2)]->Fill(phiMet[0],weight);
  for(int i_n=0; i_n<histos[int(histos.size()-1)]->GetNbinsX(); i_n++)
    if(njets>=i_n) 
      histos[int(histos.size()-1)]->AddBinContent(i_n+1,weight);  
}

void MuonStudy::Loop(int nevents) {
  if(fChain == 0) return;

  vector< vector<TH1D*> > Histos;
  vector< vector<TH2D*> > Histos2D;
  vector< vector<THnSparseD*> > HistosND;
 
  HistosND.push_back(CreateHistosND("MuonW_IC_calo"));
  HistosND.push_back(CreateHistosND("MuonWtau_IC_calo"));
  HistosND.push_back(CreateHistosND("OtherMuon_IC_calo"));
  HistosND.push_back(CreateHistosND("MuonW_kT_calo"));
  HistosND.push_back(CreateHistosND("MuonWtau_kT_calo"));
  HistosND.push_back(CreateHistosND("OtherMuon_kT_calo"));
  HistosND.push_back(CreateHistosND("MuonW_SIS_calo"));
  HistosND.push_back(CreateHistosND("MuonWtau_SIS_calo"));
  HistosND.push_back(CreateHistosND("OtherMuon_SIS_calo"));
  HistosND.push_back(CreateHistosND("MuonW_IC_track"));
  HistosND.push_back(CreateHistosND("MuonWtau_IC_track"));
  HistosND.push_back(CreateHistosND("OtherMuon_IC_track"));
  HistosND.push_back(CreateHistosND("MuonW_kT_track"));
  HistosND.push_back(CreateHistosND("MuonWtau_kT_track"));
  HistosND.push_back(CreateHistosND("OtherMuon_kT_track"));
  HistosND.push_back(CreateHistosND("MuonW_SIS_track"));
  HistosND.push_back(CreateHistosND("MuonWtau_SIS_track"));
  HistosND.push_back(CreateHistosND("OtherMuon_SIS_track"));

  Histos2D.push_back(CreateHistos2D("MuonW_IC_calo"));
  Histos2D.push_back(CreateHistos2D("MuonWtau_IC_calo"));
  Histos2D.push_back(CreateHistos2D("OtherMuon_IC_calo"));
  Histos2D.push_back(CreateHistos2D("MuonW_kT_calo"));
  Histos2D.push_back(CreateHistos2D("MuonWtau_kT_calo"));
  Histos2D.push_back(CreateHistos2D("OtherMuon_kT_calo"));
  Histos2D.push_back(CreateHistos2D("MuonW_SIS_calo"));
  Histos2D.push_back(CreateHistos2D("MuonWtau_SIS_calo"));
  Histos2D.push_back(CreateHistos2D("OtherMuon_SIS_calo"));
  Histos2D.push_back(CreateHistos2D("MuonW_IC_track"));
  Histos2D.push_back(CreateHistos2D("MuonWtau_IC_track"));
  Histos2D.push_back(CreateHistos2D("OtherMuon_IC_track"));
  Histos2D.push_back(CreateHistos2D("MuonW_kT_track"));
  Histos2D.push_back(CreateHistos2D("MuonWtau_kT_track"));
  Histos2D.push_back(CreateHistos2D("OtherMuon_kT_track"));
  Histos2D.push_back(CreateHistos2D("MuonW_SIS_track"));
  Histos2D.push_back(CreateHistos2D("MuonWtau_SIS_track"));
  Histos2D.push_back(CreateHistos2D("OtherMuon_SIS_track"));

  Histos.push_back(CreateHistos1D("MuonW_IC_calo"));
  Histos.push_back(CreateHistos1D("MuonWtau_IC_calo"));
  Histos.push_back(CreateHistos1D("OtherMuon_IC_calo"));
  Histos.push_back(CreateHistos1D("MuonW_kT_calo"));
  Histos.push_back(CreateHistos1D("MuonWtau_kT_calo"));
  Histos.push_back(CreateHistos1D("OtherMuon_kT_calo"));
  Histos.push_back(CreateHistos1D("MuonW_SIS_calo"));
  Histos.push_back(CreateHistos1D("MuonWtau_SIS_calo"));
  Histos.push_back(CreateHistos1D("OtherMuon_SIS_calo"));
  Histos.push_back(CreateHistos1D("MuonW_IC_track"));
  Histos.push_back(CreateHistos1D("MuonWtau_IC_track"));
  Histos.push_back(CreateHistos1D("OtherMuon_IC_track"));
  Histos.push_back(CreateHistos1D("MuonW_kT_track"));
  Histos.push_back(CreateHistos1D("MuonWtau_kT_track"));
  Histos.push_back(CreateHistos1D("OtherMuon_kT_track"));
  Histos.push_back(CreateHistos1D("MuonW_SIS_track"));
  Histos.push_back(CreateHistos1D("MuonWtau_SIS_track"));
  Histos.push_back(CreateHistos1D("OtherMuon_SIS_track"));
  

  // vector<TH1D*> HMuonW1D     = CreateHistos1D("MuonW");
//   vector<TH1D*> HMuonWtau1D     = CreateHistos1D("MuonWtau");
//   vector<TH1D*> HOtherMuon1D = CreateHistos1D("OtherMuon");
  
  

  
  // vector<TH2D*> HMuonW2D     = CreateHistos("MuonW");
//   vector<TH2D*> HMuonWtau2D     = CreateHistos("MuonWtau");
//   vector<TH2D*> HOtherMuon2D = CreateHistos("OtherMuon");
  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  if(nentries > nevents) nentries = nevents;
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100 == 0)
      cout << ">>> Processing event # " << jentry << endl;

    // select between Wjets Zjets or ttbarjets
    if(AlpgenIdSelection(genAlpgenID,"Wjets") == false) continue;
    weight = genWeight; // times your k factor

    // if running on ppMuX comment-out the two lines above 
    // and uncomment the line below
    // weight = 1.;

    // Look for the gen-level muon from W
    iGenmu = -99;
    for(int i=0; i<nMc; i++) 
      if(statusMc[i] == 3) // stable particle
	if(abs(idMc[i]) == 13) // it is a muon
	  if(abs(idMc[mothMc[i]]) == 24) // its mother is a W
            if(fabs(etaMc[i]) < 2.4) // Muon within the acceptance
              if(pMc[i]*sin(thetaMc[i])> 5.) // minimum pT to make it up to the muon statistions
		iGenmu = i;

    // Look for the gen-level muon from tau (from W)
    // magari domani

     iGenTaumu = -99;
     for(int i=0; i<nMc; i++)
       if(statusMc[i] == 3) // stable particle
         if(abs(idMc[i]) == 15) // it is a muon
           if(abs(idMc[mothMc[i]]) == 24) // its mother is a tau
	     if(fabs(etaMc[i]) < 2.4) // Muon within the acceptance
	       if(pMc[i]*sin(thetaMc[i])> 5.) // minimum pT to make it up to the muon statistions
		 iGenTaumu = i;


    iRecomu = -99;
    if(iGenmu != -99) {
      // Look for the corresponding reco mu
      dRMin = 99999999.;
      for(int i=0; i<nMuon; i++) {
	double dR = DeltaR(etaMuon[i],phiMuon[i],etaMc[iGenmu],phiMc[iGenmu]);
	if(dR < dRMin && dR<0.1) {
	  dRMin = dR;
	  iRecomu = i;
	}
      }
    } else if(iGenTaumu != -99){
      // Look for the corresponding reco mu
      dRMin = 99999999.;
      
      for(int i=0; i<nMuon; i++) {
        double dR = DeltaR(etaMuon[i], phiMuon[i], etaMc[iGenTaumu], phiMc[iGenTaumu]);
        if(dR < dRMin && dR<0.3) {
          dRMin = dR;
          iRecomu = i;
        }
      }
    }
    

    // RecHit thresholds for CaloTowers
    vector<float> thresh;
    thresh.push_back(.9);       //HB
    thresh.push_back(1.1);      //etc...
    thresh.push_back(1.4);   
    thresh.push_back(1.4);
    thresh.push_back(1.2);
    thresh.push_back(1.8);
    thresh.push_back(.09);
    thresh.push_back(.45);
    thresh.push_back(.2);
    thresh.push_back(.45);

    CaloF = 0.06;
    double jet_et_cut = 30;
    double jet_pt_cut = 15;
    if(nPV < 1) continue;

    iPV = -1;
    double maxpt = -1;
    for(int i = 0; i < nPV; i++){
      if(SumPtPV[i] > maxpt){
	maxpt = SumPtPV[i];
	iPV = i;
      }
    }

    double zPV = 0.0;
    if(nPV > 0) zPV = PVzPV[iPV];
    // 3 := making CaloTowers with depth from 'CaloF' depth (%)
    vector<CaloTower> calotowers = CreateCaloTowers(thresh, zPV, 3);

    int nfound = 0;
    int iCT = 0;
    while(nfound == 0 && iCT< int(calotowers.size())) {
      if(calotowers[iCT].et() > 0.5) nfound = 1;
      iCT++;
    }

    if(nfound == 1) {

      //iterative cone jets = jet_IC_calo
      vector<Jet> jet_IC_temp = SortJet(CMSIterativeConeAlgorithm(calotowers, 0.5, 1.0));
      //Fast kT jets = jet_kT_calo
      vector<Jet> jet_kT_temp = SortJet(FastJetAlgorithm(calotowers, 0.4, 1.0));
      //SISCone = jet_SIS_calo
      vector<Jet> jet_SIS_temp = SortJet(SISCone(calotowers, 0.5, 0.0));
      
      jet_SIS_track.clear();
      jet_kT_track.clear();
      jet_IC_track.clear();
      jet_SIS_calo.clear();
      jet_kT_calo.clear();
      jet_IC_calo.clear();
            
      for(int i = 0; i < jet_IC_temp.size(); i++){
	if(jet_IC_temp[i].et() > jet_et_cut && fabs(jet_IC_temp[i].eta()) <= 3.0) jet_IC_calo.push_back(jet_IC_temp[i]);
      }

      for(int i = 0; i < jet_kT_temp.size(); i++){
	if(jet_kT_temp[i].et() > jet_et_cut && fabs(jet_kT_temp[i].eta()) <= 3.0) jet_kT_calo.push_back(jet_kT_temp[i]);
      }
      
      for(int i = 0; i < jet_SIS_temp.size(); i++){
	if(jet_SIS_temp[i].et() > jet_et_cut && fabs(jet_SIS_temp[i].eta()) <= 3.0) jet_SIS_calo.push_back(jet_SIS_temp[i]);
      }

      jet_kT_temp.clear();
      jet_SIS_temp.clear();
      jet_IC_temp.clear();

      // W Muon Histograms
      if(iGenmu != -99) { // The W muon was found
	FillHistos(Histos[0], Histos2D[0], HistosND[0], iRecomu, jet_IC_calo);
	FillHistos(Histos[3], Histos2D[3], HistosND[3], iRecomu, jet_kT_calo);
	FillHistos(Histos[6], Histos2D[6], HistosND[6], iRecomu, jet_SIS_calo);
      }
      
      // W Muon from tau Histograms
      if(iGenTaumu != -99) { // The W muon was found
	FillHistos(Histos[1], Histos2D[1], HistosND[1], iRecomu, jet_IC_calo);
	FillHistos(Histos[4], Histos2D[4], HistosND[4], iRecomu, jet_kT_calo);
	FillHistos(Histos[7], Histos2D[7], HistosND[7], iRecomu, jet_SIS_calo);
      }
      
      // Other Muons Histograms
      for(int i=0;i<nMuon; i++) {
	if(i != iRecomu) { // This is not the W muon
	  FillHistos(Histos[2], Histos2D[2], HistosND[2], i, jet_IC_calo);
	  FillHistos(Histos[5], Histos2D[5], HistosND[5], i, jet_kT_calo);
	  FillHistos(Histos[8], Histos2D[8], HistosND[8], i, jet_SIS_calo);
	}
      }
    }
    
    calotowers.clear();
    
    vector<CaloTower> track_collection;
    int iTk;
    if(iRecomu > -1){
      iTk = FindTrackMu(iRecomu);
    } else {
      iTk = -1;
    }
    TVector3 vPV(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
    for(int i = 0; i < nTrack; i++){
      TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
      vT = vT-vPV;
      //remove reco muon from track collection
      if(i == iTk) continue;
      if(fabs(vT.z()) > 0.1) continue;
      if(fabs(vT.Mag()) > 1.0) continue;
      //if(abs(chargeTrack[i]) > 1) continue;
      double pt = sqrt(pxTrack[i]*pxTrack[i]+pyTrack[i]*pyTrack[i]);
      if(pt < 1.2) continue;
      if(pt > 500.) continue;
      if(trackNormalizedChi2Track[i] > 20.0) continue;
      if(fabs(trackDxyPVTrack[i]) > .6) continue;
      if(fabs(etaTrack[i]) > 2.4) continue;
      if(trackValidHitsTrack[i] < 5) continue;
      TVector3 v;
      v.SetPtEtaPhi(pt, etaTrack[i], phiTrack[i]);
      track_collection.push_back(CaloTower(pt*cosh(etaTrack[i]), 0.0, v, v, v));
    }

    nfound = 0;
    iCT = 0;
    while(nfound == 0 && iCT < int(track_collection.size())) {
      if(track_collection[iCT].et() > 0.5) nfound = 1;
      iCT++;
    }
    
    if(nfound == 1) {
     
      //IC jets (track) = jet_IC_track
      vector<Jet> jet_IC_temp = SortJet(CMSIterativeConeAlgorithm(track_collection, 0.5, 1.0));
      //Fast kT jets = jet_kT_track
      vector<Jet> jet_kT_temp = SortJet(FastJetAlgorithm(track_collection, 0.4, 1.0));
      //SISCone = jet_SIS_track
      vector<Jet> jet_SIS_temp = SortJet(SISCone(track_collection, 0.5, 0.0));
      
      for(int i = 0; i < jet_IC_temp.size(); i++){
	if(jet_IC_temp[i].pt() > jet_pt_cut) jet_IC_track.push_back(jet_IC_temp[i]);
      }
      for(int i = 0; i < jet_kT_temp.size(); i++){
	if(jet_kT_temp[i].pt() > jet_pt_cut) jet_kT_track.push_back(jet_kT_temp[i]);
      }
      for(int i = 0; i < jet_SIS_temp.size(); i++){
	if(jet_SIS_temp[i].pt() > jet_pt_cut) jet_SIS_track.push_back(jet_SIS_temp[i]);
      }
      
      jet_kT_temp.clear();
      jet_SIS_temp.clear();
      jet_IC_temp.clear();

      // fill plots as a function of jet multiplicity
      
      // W Muon Histograms
      if(iGenmu != -99) { // The W muon was found
	FillHistos(Histos[9], Histos2D[9], HistosND[9], iRecomu, jet_IC_track);
	FillHistos(Histos[12], Histos2D[12], HistosND[12], iRecomu, jet_kT_track);
	FillHistos(Histos[15], Histos2D[15], HistosND[15], iRecomu, jet_SIS_track);
      }
      
      // W Muon from tau Histograms
      if(iGenTaumu != -99) { // The W muon was found
	FillHistos(Histos[10], Histos2D[10], HistosND[10], iRecomu, jet_IC_track);
	FillHistos(Histos[13], Histos2D[13], HistosND[13], iRecomu, jet_kT_track);
	FillHistos(Histos[16], Histos2D[16], HistosND[16], iRecomu, jet_SIS_track);	
      }
      
      // Other Muons Histograms
      for(int i=0;i<nMuon; i++) {
	if(i != iRecomu) { // This is not the W muon
	  FillHistos(Histos[11], Histos2D[11], HistosND[11], i, jet_IC_track);
	  FillHistos(Histos[14], Histos2D[14], HistosND[14], i, jet_kT_track);
	  FillHistos(Histos[17], Histos2D[17], HistosND[17], i, jet_SIS_track);       	  
	}
      }

    }
    
    track_collection.clear();
    
    jet_SIS_track.clear();
    jet_kT_track.clear();
    jet_IC_track.clear();
    jet_SIS_calo.clear();
    jet_kT_calo.clear();
    jet_IC_calo.clear();
  }

  TFile *file = new TFile(outfilename.c_str(),"RECREATE");
  
  // WriteHistos(HMuonW2D, file, "MuonW2D");
  //   WriteHistos(HMuonWtau2D, file, "MuonWtau2D");
  //   WriteHistos(HMuonW1D, file, "MuonW1D");
  //   WriteHistos(HMuonWtau1D, file, "MuonWtau1D");
  
  //   WriteHistos(HOtherMuon2D, file, "OtherMuon2D");
  //   WriteHistos(HOtherMuon1D, file, "Other
  
  WriteHistos(Histos[0],  file, "MuonW_IC_calo");
  WriteHistos(Histos[1],  file, "MuonWtau_IC_calo");
  WriteHistos(Histos[2],  file, "OtherMuon_IC_calo");
  WriteHistos(Histos[3],  file, "MuonW_kT_calo");
  WriteHistos(Histos[4],  file, "MuonWtau_kT_calo");
  WriteHistos(Histos[5],  file, "OtherMuon_kT_calo");
  WriteHistos(Histos[6],  file, "MuonW_SIS_calo");
  WriteHistos(Histos[7],  file, "MuonWtau_SIS_calo");
  WriteHistos(Histos[8],  file, "OtherMuon_SIS_calo");
  WriteHistos(Histos[9],  file, "MuonW_IC_track");
  WriteHistos(Histos[10], file, "MuonWtau_IC_track");
  WriteHistos(Histos[11], file, "OtherMuon_IC_track");
  WriteHistos(Histos[12], file, "MuonW_kT_track");
  WriteHistos(Histos[13], file, "MuonWtau_kT_track");
  WriteHistos(Histos[14], file, "OtherMuon_kT_track");
  WriteHistos(Histos[15], file, "MuonW_SIS_track");
  WriteHistos(Histos[16], file, "MuonWtau_SIS_track");
  WriteHistos(Histos[17], file, "OtherMuon_SIS_track");

  WriteHistos(Histos2D[0],  file, "MuonW_IC_calo");
  WriteHistos(Histos2D[1],  file, "MuonWtau_IC_calo");
  WriteHistos(Histos2D[2],  file, "OtherMuon_IC_calo");
  WriteHistos(Histos2D[3],  file, "MuonW_kT_calo");
  WriteHistos(Histos2D[4],  file, "MuonWtau_kT_calo");
  WriteHistos(Histos2D[5],  file, "OtherMuon_kT_calo");
  WriteHistos(Histos2D[6],  file, "MuonW_SIS_calo");
  WriteHistos(Histos2D[7],  file, "MuonWtau_SIS_calo");
  WriteHistos(Histos2D[8],  file, "OtherMuon_SIS_calo");
  WriteHistos(Histos2D[9],  file, "MuonW_IC_track");
  WriteHistos(Histos2D[10], file, "MuonWtau_IC_track");
  WriteHistos(Histos2D[11], file, "OtherMuon_IC_track");
  WriteHistos(Histos2D[12], file, "MuonW_kT_track");
  WriteHistos(Histos2D[13], file, "MuonWtau_kT_track");
  WriteHistos(Histos2D[14], file, "OtherMuon_kT_track");
  WriteHistos(Histos2D[15], file, "MuonW_SIS_track");
  WriteHistos(Histos2D[16], file, "MuonWtau_SIS_track");
  WriteHistos(Histos2D[17], file, "OtherMuon_SIS_track");

  WriteHistosND(HistosND[0],  file, "MuonW_IC_calo");
  WriteHistosND(HistosND[1],  file, "MuonWtau_IC_calo");
  WriteHistosND(HistosND[2],  file, "OtherMuon_IC_calo");
  WriteHistosND(HistosND[3],  file, "MuonW_kT_calo");
  WriteHistosND(HistosND[4],  file, "MuonWtau_kT_calo");
  WriteHistosND(HistosND[5],  file, "OtherMuon_kT_calo");
  WriteHistosND(HistosND[6],  file, "MuonW_SIS_calo");
  WriteHistosND(HistosND[7],  file, "MuonWtau_SIS_calo");
  WriteHistosND(HistosND[8],  file, "OtherMuon_SIS_calo");
  WriteHistosND(HistosND[9],  file, "MuonW_IC_track");
  WriteHistosND(HistosND[10], file, "MuonWtau_IC_track");
  WriteHistosND(HistosND[11], file, "OtherMuon_IC_track");
  WriteHistosND(HistosND[12], file, "MuonW_kT_track");
  WriteHistosND(HistosND[13], file, "MuonWtau_kT_track");
  WriteHistosND(HistosND[14], file, "OtherMuon_kT_track");
  WriteHistosND(HistosND[15], file, "MuonW_SIS_track");
  WriteHistosND(HistosND[16], file, "MuonWtau_SIS_track");
  WriteHistosND(HistosND[17], file, "OtherMuon_SIS_track");

}

double MuonStudy::pT(double px, double py) {
  return sqrt(px*px+py*py);
}

int MuonStudy::FindTrackMu(int MuInd) {
  int ind = -99;
  double dRmin_jet = 99999999999.;
  for(int i=0; i< nTrack; i++) {
    double dR = DeltaR(etaMuon[MuInd],phiMuon[MuInd],etaTrack[i],phiTrack[i]);
    if(dR<dRmin_jet) {
      ind = i;
      dRmin_jet = dR;
    }
  }

  if(dRmin_jet > 0.01) ind = -1;
  
  //  if(ind == -99) {
  //    cout << " I did not find the track" << endl;
  //    cout << dRmin << endl;
  //  }

  return ind;
}

double MuonStudy::mT(int imuon) {
  double PtVectProduct=
    pxMuon[imuon]*(pxMet[0]-pxMuon[imuon])+
    pyMuon[imuon]*(pyMet[0]-pyMuon[imuon]);
  double cosphi= PtVectProduct/pT(pxMuon[imuon],pyMuon[imuon])/pT(pxMet[0]-pxMuon[imuon],pyMet[0]-pyMuon[imuon]);
  return sqrt(2.*pT(pxMuon[imuon],pyMuon[imuon])*pT(pxMet[0]-pxMuon[imuon],pyMet[0]-pyMuon[imuon])*(1-cosphi));
}
