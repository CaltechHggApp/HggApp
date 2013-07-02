//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

/// The MuonStudy class can be used to perform fast check
/// on input ntuples (in a format compatible to VecbosBase)

#ifndef MuonStudy_h
#define MuonStudy_h

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <vector>
#include <THnSparse.h>

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"

using namespace std;

class MuonStudy : public Vecbos{
public:

  MuonStudy(TTree *tree=0, string outname="test.root"); /// Class Constructor
  virtual ~MuonStudy();     /// Class Destructor
  /// The function to run on each events
  void Loop(int);
    
private:

  string outfilename;

  double mT(int imuon);

  /// Create a set of histograms. Takes as input the name of the 
  /// directory were to write the histograms, also used to give 
  /// names to the histograms (to avoid memory problems if used
  /// for more than a set of histograms)
  vector<TH2D*> CreateHistos2D(string dirname); 
  vector<TH1D*> CreateHistos1D(string dirname); 
  vector<THnSparseD*> CreateHistosND(string dirname);
  void WriteHistosND(vector<THnSparseD*> histos, TFile* file, string dirname);

  /// Finf the track corresponding to the muon
  int FindTrackMu(int MuInd);

  double pT(double px, double py);

  double weight;

  /// Fill the histograms passes as input
  void FillHistos(vector<TH1D*> histos,vector<TH2D*> histos2d, vector<THnSparseD*> histosNd, int iMu, vector<Jet> jets );
  
  double dRMin;
  int iRecomu, iGenmu, iGenTaumu, njets;
  vector<Jet> jet_IC_calo, jet_SIS_calo, jet_kT_calo;
  vector<Jet> jet_IC_track, jet_SIS_track, jet_kT_track;
};
#endif
