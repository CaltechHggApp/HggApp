//--------------------------------------------------------------
// Description:
//    Auxiliary class for selection of reconstructed Vecbos+jets
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//    Maurizio Pierini 
//    CERN
//    Thiago Tomei 
//    SPRACE, Sao Paulo, Brazil
//    Chris Rogan
//    Caltech
//    Yi Chen (adding electron/muon function, b-tag working points and hcal noise wrapper function)
//--------------------------------------------------------------

/// The Vecbos class is an auxiliary class which contains basic
/// functionality useful for any analysis of Vecbos+jets events.
/// It derives from VecbosBase.

#ifndef Vecbos_h
#define Vecbos_h

// uncomment if tree content includes the calotowers block
// in order to use methods using them
// #define USECALOTOWERS

// std includes
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include <numeric>
#include <list>
#include "combination.hh"

// FASTJET includes
#include "FASTJET/include/fastjet/PseudoJet.hh"
#include "FASTJET/include/fastjet/ClusterSequence.hh"
#include "FASTJET/include/fastjet/ClusterSequenceActiveArea.hh"
#include "FASTJET/include/fastjet/SISConePlugin.hh"

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TVector.h>
#include <TLorentzVector.h>

#include "CommonTools/include/Utils.hh"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"

// VecbosApp includes
#include "Jet.hh"
#include "BTagJet.hh"
#include "MET.hh"
#include "CaloTower.hh"
#include "CoolTools.hh"
#include "VecbosBase.hh"
#include "JetCorrectionUncertainty.h"

using namespace std;

/// The Vecbos class is an auxiliary class which contains basic
/// functionality useful for any analysis of Vecbos+jets events.
/// It derives from VecbosBase.
/// More specialized analysis classes should derive from Vecbos.
class Vecbos : public VecbosBase
{
public:
  typedef std::pair<unsigned int,unsigned int> aLSSegment;
  typedef std::vector< std::pair<unsigned int,unsigned int> > LSSegments;
  typedef unsigned int aRun;
  typedef std::map< aRun, LSSegments > runsLSSegmentsMap;
  typedef std::pair < aRun, LSSegments > aRunsLSSegmentsMapElement;

  /// This function is needed to read the setting parameters
  /// from an input text file.
  /// Reading line by line the input text file, create a
  /// map that relates each variable with its value
  void ReadParameters(const char* );

  /// Set verbosity level.
  void SetVerbose(bool v){verbose = v;}
  /// Class Constructor
  Vecbos(TTree *tree=0);
  /// Class Destructor
  virtual ~Vecbos();

  /// The primary vertex of the event.
  Int_t iPV;
  /// The depth tor the shower in the calorimeter
  Double_t CaloF;

  /// Fill RunLSMap according to json file
  void fillRunLSMap();
  /// Set Good Run LS
  void setJsonGoodRunList(const string& jsonFilePath);
  /// check if Run/LS is a good one
  bool isGoodRunLS();
  /// reload TriggerMask if necessary (data file is changed). Should be called for each event inside the event loop
  bool reloadTriggerMask(bool newVersion=false);
  bool reloadTriggerMask(int runN);
  /// Gives the right trigger for the desired run range
  std::string getHLTPathForRun(int runN, std::string fullname);
  /// set the list of required trigger to produce the bitmask
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers);
  //check if the event passed HLT. To be called per event
  bool hasPassedHLT();
  //check for matching HLT object
  bool triggerMatch(float eta, float phi, float Dr);
  //get the value of the requested bits
  vector<int> getHLTOutput();

private:
  struct JetConfig;
  JetConfig *theJetConfig;
  
  ///goodRUN/LS list
  runsLSSegmentsMap goodRunLS; 
  std::string jsonFile;

  std::string lastFile;
  std::vector<std::string> requiredTriggers;

  JetCorrectionUncertainty *jecUnc_calo;
  JetCorrectionUncertainty *jecUnc_PF;

protected:

  // compute M_R
  double CalcMR(TLorentzVector ja, TLorentzVector jb);
  // compute MR*
  double CalcMRstar(TLorentzVector ja, TLorentzVector jb);
  // compute gamma*MR*
  double CalcGammaMRstar(TLorentzVector ja, TLorentzVector jb);
  // compute M_R'
  double CalcMRP(TLorentzVector ja, TLorentzVector jb, TVector3 met); 
  // compute M_T^R
  double CalcMTR(TLorentzVector ja, TLorentzVector jb, TVector3 met);

  //! the list of required triggers
  std::vector<int> m_requiredTriggers;
  
  bool AlpgenIdSelection(double alpgenid, string sample); ///< Filter events according to the ALPGEN ID

  bool verbose; ///< Verbosity of printouts.

  /// Write the given histograms in the output file, storing them
  /// in a directory (named as specified)
  void WriteHistos(vector<TH1D*> histos, TFile* file, string dirname);
  void WriteHistos(vector<TProfile*> histos, TFile* file, string dirname);
  void WriteHistos(vector<TH2D*> histos, TFile* file, string dirname);

  /// Taking as input a map containing the name and the
  /// variables of the setting inputs, it initializes the
  /// analysis parameters. The analysis parameters
  /// are private members of the VecbosBase class.
  virtual void AssignParameters(map<string, double>);

  /// Like AssignParameters(map<string, double>)
  /// but for initialization of string parameters
  virtual void AssignParameters(map<string, string>);
 
  /// Initialize the analysis parameters to their
  /// default values. All the values can be changed by
  /// the user specifying their value in a config file,
  /// to pass as input to the ReadParameters() function.
  virtual void InitParameters();

  // Analysis Parameters

  /// String identifying the process to run on
  /// Set to one between: "Wjets", "Zjets", "ttbarjets"
  /// Default value is "Zjets"
  string process;

  /// Maximum rapidity value, defining the barrel region.
  /// Default value is 1.3
  double barrellimit; 

  /// Maximum rapidity value, defining the endcap region.
  /// The Endcap region goes from barrellimit to endcaplimit.
  /// Default value is 3.0
  double endcaplimit; 

  // Some useful templates. The comparators work both with Jets and
  // with TLorentzVectors. Some work with TVector2, TVector3 as well.
  template <typename T> 
  T DeltaPhi(T phi1, T phi2); //< Delta phi in radians in between two angles.
  template <typename T>
  T DeltaR(T eta1, T phi1, T eta2, T phi2); //< Delta R in between two pairs (eta,phi).

  /// Creates a set of TLorentzVectors from the particles in the event, with a given status.
  vector<TLorentzVector> ParticlesFromMc(int status); 
  /// Creates a set of TLorentzVector from the particles in the event, with a given status
  /// and with a given set of IDs. Useful to get, e.g., all muons in the event.
  /// You HAVE to specify both positive and negative IDs.
  vector<TLorentzVector> ParticlesFromMcWithId(int status, const vector<int>& allowed); 
  /// Creates a set of TLorentzVector from the particles in the event, with a given status
  /// and NOT with a given set of IDs. Useful to get, e.g., all particles but muons in the event.
  /// You HAVE to specify both positive and negative IDs.
  vector<TLorentzVector> ParticlesFromMcWithNotId(int status, const vector<int>& forbidden);
  /// Creates a set of TLorentzVector from the stable charged particles in the event.
  vector<TLorentzVector> ParticlesFromMcCharged();
  /// Creates a set of TLorentzVector from the tracks - with an optional cut in pt.
  vector<TLorentzVector> Tracks(double thePtCut = 0.);
  
  /// Creates a set of TLorentzVector close to another TLorentzVector by deltaR
  vector<TLorentzVector> CloseInEtaPhi(const vector<TLorentzVector>& set, const TLorentzVector& v, double maxDistance);

#ifdef USECALOTOWERS
  // Internal functions for jet reclustering in the Vecbos code.
  // If you use these functions, then you should use the jets you 
  // create instead of the ones in the ntuples.
  /// Standard set of calorimetric thresholds.
  virtual vector<float> DefaultCaloThresholds();
  /// Creates a set of CaloTowers.
  vector<CaloTower> CreateCaloTowers(vector<float>, float, int);
  /// Create MET from the CaloTowers.
  MET CreateMET(vector<CaloTower>);
  /// Create MET from Jets.
  MET CreateMET(vector<Jet>);
  /// Create MET from TLorentzVector
  MET CreateMET(vector<TLorentzVector>);
#endif

  /// Creates a set of IC jets from CaloTowers.
  vector<Jet> CMSIterativeConeAlgorithm(vector<CaloTower>,double,double);
  /// Creates a set of IC jets from 4vectors. 
  vector<Jet> CMSIterativeConeAlgorithm(vector<TLorentzVector>,double,double); 
  /// Creates a set of KT jets from CaloTowers.
  vector<Jet> FastJetAlgorithm(vector<CaloTower>,double,double);
  /// Creates a set of KT jets from 4vectors.
  vector<Jet> FastJetAlgorithm(vector<TLorentzVector>,double,double);
  /// Creates a set of SISCone jets from CaloTowers.
  vector<Jet> SISCone(vector<CaloTower>,double,double);
  /// Creates a set of SISCone jets from 4vectors.
  vector<Jet> SISCone(vector<TLorentzVector>,double,double);

  // Useful jet functions.
  /// Creates a vector<pair<Jet,Jet>, for matching jets.
  vector<pair<Jet,Jet> > OneToOneMatch(vector<Jet>,vector<Jet>);
  /// Sorts a vector of Jets by Pt, although the function name doesn't say so... 
  vector<Jet> SortJet(vector<Jet> v);
  /// Sorts a vector of Jets by Et.
  vector<Jet> SortJetByEt(vector<Jet> v);
  /// Correct a jet, given a vertex.
  Jet CorrectJet(Jet J, double vtx);
  /// Give the event EM fraction, given a collection of jets and cuts.
  double EventEMF(vector<Jet> vJ, double ptCut, double etaCut);
  /// Give the jet charged fraction.
  double JetCHF(Jet, double, int&);
  /// Give the event charged fraction, given a collection of jets and cuts.
  double EventCHF(vector<Jet> vJ, double ptCut, double etaCut);
  /// Read the gen jets from the event.
  vector<Jet> GetGenJets();
  /// Read the uncorrected jets from the event.
  vector<Jet> GetUncorrJets();
  /// Read the corrected jets from the event.
  /// scale energy = apply JES uncertainty: +/-1 * uncertainty; 0 = do not apply (default)
  vector<Jet> GetCorrJets(int scaleEnergy=0);
  /// Read the btag algorithms output
  vector<BTagJet> GetBTagCorrJets();
  /// Read the uncorrected PF jets from the event.
  vector<Jet> GetUncorrPFJets();
  /// Read the corrected PF jets from the event. 
  /// scale energy = apply JES uncertainty: +/-1 * uncertainty; 0 = do not apply (default)
  vector<Jet> GetCorrPFJets(int scaleEnergy=0);
  /// Read the PF btag algorithms output
  vector<BTagJet> GetBTagCorrPFJets();
  /// Read the uncorrected PF jets from the event (with Fastjet PU,UE subtraction).
  vector<Jet> GetUncorrPUPFJets();
  /// Read the corrected PF jets from the event (with Fastjet PU,UE subtraction). 
  /// scale energy = apply JES uncertainty: +/-1 * uncertainty; 0 = do not apply (default)
  vector<Jet> GetCorrPUPFJets(int scaleEnergy=0);
  /// Read the PF btag algorithms output (with Fastjet PU,UE subtraction)
  vector<BTagJet> GetBTagCorrPUPFJets();

  /// Get pt given x/y coordinates
  float GetPt(float px, float py) { return TMath::Sqrt(px*px + py*py); }
  /// Combinatoric matching stuff
  std::vector<int> AlgoBruteForce(int, int);
  std::vector<int> AlgoSwitchMethod(int, int);
  double lengtht( std::vector<int> );
  std::vector< std::vector<float> > AllDist;

  // useful electron functions
  /// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
  float SigmaiEiE(int electron);
  /// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
  float SigmaiPiP(int electron);
  // get the likelihood electron ID
  float likelihoodRatio(int eleIndex, ElectronLikelihood &lh);
  /// bremsstrahlung fraction
  float FBrem(int electron);
  /// e9Esc
  float E9ESC(int electron);
  /// e25Esc
  float E25ESC(int electron);
  /// esc
  float ESC(int electron);
  /// eseed
  float Eseed(int electron);
  /// dxy and dsz parameters with respect to PV for electron tracks
  double eleDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);
  double eleDszPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);
  /// evaluate di-jet mass (only in events with =2 jets, return -1 otherwise)
  double GetDiJetMass(std::vector<Jet> jets);
  /// evaluate di-jet pT (only in events with =2 jets, return -1 otherwise)
  double GetDiJetPt(std::vector<Jet> jets);
  /// evaluate delta-eta between the two jets
  double GetDeltaEtaJet(std::vector<Jet> jets);
  /// evaluate delta-phi between jet and met
  double GetDeltaPhiJetMet(std::vector<Jet> jets, TVector3 met);

  /// returns the pt of the true photon in a photon + jet event (approx pthat)
  /// needed to make exclusive photon+jets samples (retain pT_sample(i) < pthat < pT_sample(i+1) )
  float photonPt();

  // Check if electron passes id working point 95/90/85/80/70/60
  // Following https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID2011 (2011 May 25)
  bool electronPassWP(int index, int WP);
  bool electronPassWP95(int index);
  bool electronPassWP90(int index);
  bool electronPassWP85(int index);
  bool electronPassWP80(int index);
  bool electronPassWP70(int index);
  bool electronPassWP60(int index);

  // taken from the VBTF_new.* from Chris
  bool muonPassLoose(int index);
  bool muonPassTight(int index);

  // Hcal noise bit wrapper
  bool eventPassHcalFilter();

  // b-tag working points for jets
  bool caloJetPassTCHEL(int index);
  bool caloJetPassTCHEM(int index);
  bool caloJetPassTCHET(int index);
  bool caloJetPassTCHPL(int index);
  bool caloJetPassTCHPM(int index);
  bool caloJetPassTCHPT(int index);
  bool caloJetPassSSVHEM(int index);
  bool caloJetPassSSVHET(int index);
  bool caloJetPassSSVHPT(int index);
  
  bool pfJetPassTCHEL(int index);
  bool pfJetPassTCHEM(int index);
  bool pfJetPassTCHET(int index);
  bool pfJetPassTCHPL(int index);
  bool pfJetPassTCHPM(int index);
  bool pfJetPassTCHPT(int index);
  bool pfJetPassSSVHEM(int index);
  bool pfJetPassSSVHET(int index);
  bool pfJetPassSSVHPT(int index);
};

/// Auxiliary classes for jet reclustering.
class CaloPointMy {
public:
  CaloPointMy (double fZ, double fEta, double F) {
    // Typical fractional depth of clusters in calorimeters.
    //      const double F = 0.15;
    //const double F = 0.1;
    //const double F = CaloF;
    const double F1 = 1. - F;
    const double R_BARREL = F1*143.+F*407.; // 1/2(EBrin+HOrout) from CaloTowerHardcodeGeometryLoader
    const double Z_ENDCAP = F1*320.+F*568.; // 1/2(EEz+HEz)
    const double R_FORWARD = Z_ENDCAP / sqrt (cosh(3.)*cosh(3.0) -1.); // eta=3
    const double Z_FORWARD = 1100.+F*165.;
    const double ETA_MAX = 5.2;
    const double Z_BIG = 1.e5;
    
    if (fZ > Z_ENDCAP) fZ = Z_ENDCAP-1.;
    if (fZ < -Z_ENDCAP) fZ = -Z_ENDCAP+1; // sanity check
    
    double tanThetaAbs = sqrt (cosh(fEta)*cosh(fEta) - 1.);
    double tanTheta = fEta >= 0 ? tanThetaAbs : -tanThetaAbs;
    
    double rEndcap = tanTheta == 0 ? 1.e10 : 
      fEta > 0 ? (Z_ENDCAP - fZ) / tanTheta : (-Z_ENDCAP - fZ) / tanTheta;
    if (rEndcap > R_BARREL) { // barrel
      mR = R_BARREL;
      mZ = fZ + R_BARREL * tanTheta; 
    }
    else {
      double zRef = Z_BIG; // very forward;
      if (rEndcap > R_FORWARD) zRef = Z_ENDCAP; // endcap
      else if (fabs (fEta) < ETA_MAX) zRef = Z_FORWARD; // forward
      
      mZ = fEta > 0 ? zRef : -zRef;
      mR = fabs ((mZ - fZ) / tanTheta);
    }
  }
  
  double etaReference (double fZ) {
    TVector3  p(r(), 0., z() - fZ);
    return p.Eta();
  }
  
  double thetaReference (double fZ) {
    TVector3 p(r(), 0., z() - fZ);
    return p.Theta();
  }
  
  double z() const {return mZ;}
  double r() const {return mR;}
  
private:
  double mZ;
  double mR;
};

class CaloPointECAL {
public:
  CaloPointECAL (double fZ, double fEta) {
    // Typical fractional depth of clusters in calorimeters.
    //      const double F = 0.15;
    const double F = 0.1;
    const double F1 = 1. - F;
    const double R_BARREL = 129.0;
    const double Z_ENDCAP = F1*320.+F*568.; // 1/2(EEz+HEz)
    const double R_FORWARD = Z_ENDCAP / sqrt (cosh(3.)*cosh(3.0) -1.); // eta=3
    const double Z_FORWARD = 1100.+F*165.;
    const double ETA_MAX = 5.2;
    const double Z_BIG = 1.e5;
    
  

    if (fZ > Z_ENDCAP) fZ = Z_ENDCAP-1.;
    if (fZ < -Z_ENDCAP) fZ = -Z_ENDCAP+1; // sanity check
    
    double tanThetaAbs = sqrt (cosh(fEta)*cosh(fEta) - 1.);
    double tanTheta = fEta >= 0 ? tanThetaAbs : -tanThetaAbs;
    
    double rEndcap = tanTheta == 0 ? 1.e10 : 
      fEta > 0 ? (Z_ENDCAP - fZ) / tanTheta : (-Z_ENDCAP - fZ) / tanTheta;
    if (rEndcap > R_BARREL) { // barrel
      mR = R_BARREL;
      mZ = fZ + R_BARREL * tanTheta; 
    }
    else {
      double zRef = Z_BIG; // very forward;
      if (rEndcap > R_FORWARD) zRef = Z_ENDCAP; // endcap
      else if (fabs (fEta) < ETA_MAX) zRef = Z_FORWARD; // forward
      
      mZ = fEta > 0 ? zRef : -zRef;
      mR = fabs ((mZ - fZ) / tanTheta);
    }
  }
  
  double etaReference (double fZ) {
    TVector3  p(r(), 0., z() - fZ);
    return p.Eta();
  }
  
  double thetaReference (double fZ) {
    TVector3 p(r(), 0., z() - fZ);
    return p.Theta();
  }
  
  double z() const {return mZ;}
  double r() const {return mR;}
  
private:
  double mZ;
  double mR;
};

/// Template function has to be moved here
/// as to avoid linker errors. See rationale
/// in C++ reference
/// http://www.parashift.com/c++-faq-lite/templates.html#faq-35.12

template <typename T>
T Vecbos::DeltaPhi(T phi1, T phi2) {
  T result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

template <typename T>
T Vecbos::DeltaR(T eta1, T phi1, T eta2, T phi2) {
  T dphi = DeltaPhi(phi1,phi2);
  T result = sqrt((eta1-eta2)*(eta1-eta2)+dphi*dphi);
  return result;
}
  

#endif
