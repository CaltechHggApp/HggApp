//---------------------------------------------------------------------------
// Description:
//    Class for Zmumu candle analysis
// Authors:
//    Y. Chen
//---------------------------------------------------------------------------
#ifndef CreateWJetDataset_h
#define CreateWJetDataset_h
//---------------------------------------------------------------------------
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/TriggerMask.hh"
//---------------------------------------------------------------------------
#define DefaultCaloThreshold 30
#define DefaultUncorrectedCaloThreshold 20
#define DefaultPFThreshold 30
#define DefaultTrackThreshold 20
//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
class CreateWJetDataset;
class MuonCandidate;
struct OutputRecord;
struct JetRecord;
struct GenBHadron;
struct SingleJet;
//---------------------------------------------------------------------------
class CreateWJetDataset : public Vecbos
{
public:
  CreateWJetDataset(TTree *tree = 0); /// Class Constructor
  CreateWJetDataset(TTree *tree = 0, bool goodRunLS = false, bool isData = false,
		    int selectZMuMu = -1); /// Class Constructor
  virtual ~CreateWJetDataset();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  void SetConditions(TTree* treeCond);

private:
  double SumPt(int iMu, int iZ);
  int BestPV();
  double pTMuon(int imu);
  double DeltaPhi_PiHalf(double phi1, double phi2);

  double weight;
  vector<int> m_requiredTriggers;

  bool _isData;
  bool _goodRunLS;
  int _selectZMuMu;   // if signal MC, go into MC information and find Z and its daughters
  // -1: don't care, 1: require one, 0: no ZMuMu allowed
  TTree* _treeCond;

private:
  vector<MuonCandidate> MakeMuonCandidates();
  int NumberPassFirstLeg(vector<MuonCandidate> &Muons);
  int NumberPassSecondLeg(vector<MuonCandidate> &Muons);
  JetRecord FindCaloJets(vector<GenBHadron> & bcHadrons,vector<GenBHadron> & bcHadrons_alt,vector<GenBHadron> & bcHadrons_alt2,double Threshold = DefaultCaloThreshold, double OverallScale = 1.0, double EtaSlope = 0);
  //JetRecord FindMatchedCaloJets(vector<GenBHadron> BHadrons, double Threshold = DefaultCaloThreshold, double OverallScale = 1.0, double EtaSlope = 0);
  //JetRecord FindUncorrectedCaloJets(double Threshold = DefaultUncorrectedCaloThreshold,
  //double OverallScale = 1.0, double EtaSlope = 0);
  JetRecord FindPFJets(MuonCandidate &TheMuon, vector<GenBHadron> & bcHadrons,double Threshold = DefaultPFThreshold,
		       double OverallScale = 1.0, double EtaSlope = 0);
  //JetRecord FindTrackJets(MuonCandidate &TheMuon, int GoodPV, double Threshold = DefaultTrackThreshold,
  //		  double OverallScale = 1.0, double EtaSlope = 0);
  vector<GenBHadron> FindGenBHadrons_Alt();
  vector<GenBHadron> FindGenBHadrons_Alt2();

private:
  void MakeBranches(TTree *tree, OutputRecord *Record);

private:
  map<string, TH1D *> GenerateQM1DHistograms();
  map<string, TH2D *> GenerateQM2DHistograms();
  void WriteQMHistograms(map<string, TH1D *> Histograms);
  void WriteQMHistograms(map<string, TH2D *> Histograms);
  void DeleteQMHistograms(map<string, TH1D *> Histograms);
  void DeleteQMHistograms(map<string, TH2D *> Histograms);
  bool ExistBInGenParticle(vector<GenBHadron> & bcHadrons);
};
//---------------------------------------------------------------------------
class MuonCandidate
{
public:
  bool IsGlobal;
  bool IsPromptTight;
  bool IsTracker;
  // int PixelHit;
  // int StripHit;
  // double Chi2;
  // int ValidMuonHit;
  // int MuonStations;
  double Dxy;
  double Dz;
  double EcalIsolation;
  double HcalIsolation;
  double TrackIsolation;
  double Phi;
  double Eta;
  double PT;
  double Px;
  double Py;
  double Pz;
  int Charge;
public:
  bool PassMuonID;
  bool PassFirstLeg;
  bool PassSecondLeg;
  bool PassIsolation;
public:
  int MuonIndex;

public:
  MuonCandidate &operator =(MuonCandidate &other)
  {
    IsGlobal = other.IsGlobal;
    IsPromptTight = other.IsPromptTight;
    IsTracker = other.IsTracker;
    // PixelHit = other.PixelHit;
    // StripHit = other.StripHit;
    // Chi2 = other.Chi2;
    // ValidMuonHit = other.ValidMuonHit;
    // MuonStations = other.MuonStations;
    Dxy = other.Dxy;
    Dz = other.Dz;
    EcalIsolation = other.EcalIsolation;
    HcalIsolation = other.HcalIsolation;
    TrackIsolation = other.TrackIsolation;
    Phi = other.Phi;
    Eta = other.Eta;
    PT = other.PT;
    Px = other.Px;
    Py = other.Py;
    Pz = other.Pz;
    Charge = other.Charge;
    PassMuonID = other.PassMuonID;
    PassFirstLeg = other.PassFirstLeg;
    PassSecondLeg = other.PassSecondLeg;
    PassIsolation = other.PassIsolation;
    MuonIndex = other.MuonIndex;

    return *this;
  }
};
//---------------------------------------------------------------------------
struct JetRecord
{
public:
  int Count;
  double Eta[20];
  double Phi[20];
  double PT[20];
  double Energy[20];
  double CSV[20];
  double CSVMVA[20];
  double JBP[20];
  double JP[20];
  double TCHE[20];
  double TCHP[20];
  double SSV[20];
  double B_DR[20];
  double C_DR[20];
  double B_DR_alt[20];
  double C_DR_alt[20];
  double B_DR_alt2[20];
  double C_DR_alt2[20];

  double TotalPT;
  int CountSurvey[101];   // threshold 0, 1, ...., 100; ie., index = threshold
  int BTaggedCount;   // b-tag variable > 0.4
  int LooseBTaggedCount;   // b-tag variable > 0.7

public:
  JetRecord()
  {
    Count = 0;
    for(int i = 0; i < 20; i++)
      {
	Eta[i] = 0;
	Phi[i] = 0;
	PT[i] = 0;
	Energy[i] = 0;
	CSV[i] = 0;
	CSVMVA[i] = 0;
	JBP[i] = 0;
	JP[i] = 0;
	TCHE[i] = 0;
	TCHP[i] = 0;
	SSV[i] = 0;
	B_DR[i] = 0;
	C_DR[i] = 0;
	B_DR_alt[i] = 0;
	C_DR_alt[i] = 0;
	B_DR_alt2[i] = 0;
	C_DR_alt2[i] = 0;
      }
    TotalPT = 0;
    for(int i = 0; i < 101; i++)
      CountSurvey[i] = 0;
  }

  JetRecord &operator =(JetRecord &other)
  {
    Count = other.Count;
    for(int i = 0; i < 20; i++)
      {
	Eta[i] = other.Eta[i];
	Phi[i] = other.Phi[i];
	PT[i] = other.PT[i];
	Energy[i] = other.Energy[i];
	CSV[i] = other.CSV[i];
	CSVMVA[i] = other.CSVMVA[i];
	JBP[i] = other.JBP[i];
	JP[i] = other.JP[i];
	TCHE[i] = other.TCHE[i];
	TCHP[i] = other.TCHP[i];
	SSV[i] = other.SSV[i];
	B_DR[i] = other.B_DR[i];
	C_DR[i] = other.C_DR[i];
	B_DR_alt[i] = other.B_DR_alt[i];
	C_DR_alt[i] = other.C_DR_alt[i];
	B_DR_alt2[i] = other.B_DR_alt2[i];
	C_DR_alt2[i] = other.C_DR_alt2[i];
      }
    TotalPT = other.TotalPT;
    for(int i = 0; i < 101; i++)
      CountSurvey[i] = other.CountSurvey[i];
    BTaggedCount = other.BTaggedCount;
    LooseBTaggedCount = other.LooseBTaggedCount;

    return *this;
  }
};
//---------------------------------------------------------------------------
struct OutputRecord
{
  double RunNumber;
  double EventNumber;
  double BunchCrossing;
  double LumiSection;
  double Orbit;

  bool PassHLT;

  MuonCandidate TheMuon;
  double EcalRelativeIsolation;
  double EcalIsolationOverMT;
  double HcalRelativeIsolation;
  double HcalIsolationOverMT;
  double TrackRelativeIsolation;
  double TrackIsolationOverMT;

  JetRecord CaloJets;
  JetRecord MatchedCaloJets;
  JetRecord UncorrectedCaloJets;
  JetRecord PFJets;
  JetRecord TrackJets;

  double HighestBTag;
  double SecondHighestBTag;

  double ChrisKinematicTopVariable;   // ??????
  double MR;   // Variables designed by Chris
  double MRPrime;
  double R;
  double RPrime;

  double CaloMET;
  double PFMET;
  double TCMET;

  double CaloMT;
  double PFMT;
  double TCMT;

  bool HasB;

  double Weight;   // MC weight.  For data it's always 1
};
//---------------------------------------------------------------------------
struct GenBHadron
{
  bool Final;
  double Eta;
  double Phi;
  double Energy;
  int Status;
  int PDGID;
};
//---------------------------------------------------------------------------
#endif

