//---------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <map>

using namespace std;
//---------------------------------------------------------------------------
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
//---------------------------------------------------------------------------
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "CreateWJetDataset.hh"
//---------------------------------------------------------------------------
void CreateWJetDataset::SetConditions(TTree* treeCond)
{
  _treeCond = treeCond;
}
//---------------------------------------------------------------------------
CreateWJetDataset::CreateWJetDataset(TTree *tree, bool goodRunLS, bool isData, int selectZMuMu) : Vecbos(tree)
{
  cout << "Staring CreateWJetDataset" << endl;
  _goodRunLS = goodRunLS;
  _isData = isData;
  _selectZMuMu = selectZMuMu;

  // To read good run list!
  if(_goodRunLS && _isData)
    {
      std::string goodRunGiasoneFile = "config/vecbos/json/goodRunLS.json";
      setJsonGoodRunList(goodRunGiasoneFile);
      fillRunLSMap();
    }
}
//---------------------------------------------------------------------------
CreateWJetDataset::~CreateWJetDataset()
{
}
//---------------------------------------------------------------------------
void CreateWJetDataset::Loop(string outFileName, int start, int stop)
{
  if(fChain == NULL)
    return;

  if(stop <= start)
    return;

  if(fChain->GetEntries() == 0)
    return;

  TFile *file = new TFile(outFileName.c_str(), "RECREATE");

  // count how many events we have
  TH1D *HEventCount = new TH1D("HEventCount", "HEventCount", 1, -0.5, 0.5);

  // output record!
  OutputRecord Messenger;

  // prepare the output tree
  TTree* outTree = new TTree("WJetTree", "WJetTree");
  MakeBranches(outTree, &Messenger);
   
  // prepare quality monitoring histograms
  map<string, TH1D *> QualityMonitoring1DHistograms = GenerateQM1DHistograms();
  map<string, TH2D *> QualityMonitoring2DHistograms = GenerateQM2DHistograms();

  // start looping
  unsigned int lastLumi = 0;
  unsigned int lastRun = 0;

  Long64_t nbytes = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;


  for(Long64_t jentry = start; jentry < stop; jentry++)
    //for(Long64_t jentry = 0; jentry < 100; jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if(ientry < 0)
	break;

      if(jentry % 10000 == 0)
	cout << ">>> Processing event #" << jentry << "!!?" << endl;

      HEventCount->Fill(0);   // count total number of events

      Long64_t nb = fChain->GetEntry(jentry);
      nbytes += nb;

      // get HLT bits for the event
      /*
      if(_isData == true)
	reloadTriggerMask(true);
      else
	reloadTriggerMask(false);
      */
      reloadTriggerMask(true);
      Messenger.PassHLT = hasPassedHLT();
      
      // Good Run selection - for data
      if(_isData && _goodRunLS && !isGoodRunLS())
	{
	  if(lastRun != runNumber || lastLumi != lumiBlock)
	    {
	      lastRun = runNumber;
	      lastLumi = lumiBlock;
	      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected!" << std::endl;
	    }
	  continue;
	}

      // HLT requirement
      if(Messenger.PassHLT == false)
	continue;
      
      // At least one reconstructed muon
      if(nMuon < 1)
	continue;

      // Make muon candidates
      vector<MuonCandidate> Muons = MakeMuonCandidates();

      // We need one tight muon, and nothing else
      if(NumberPassFirstLeg(Muons) != 1)
	continue;

      // Find that muon and remove it from the muon candidate list
      MuonCandidate TheMuon;
      for(int i = 0; i < (int)Muons.size(); i++)
	{
	  if(Muons[i].PassFirstLeg == false)
            continue;
	  TheMuon = Muons[i];
	  Muons.erase(Muons.begin() + i);
	  break;
	}

      // No second muon from the rest
      if(NumberPassSecondLeg(Muons) != 0)
	continue;

      // Find the best PV
      int GoodPV = BestPV();

      // Find generator-level b hadrons
      vector<GenBHadron> bcHadrons;

      vector<GenBHadron> bcHadrons_alt;
      if(_isData == false)
	bcHadrons_alt = FindGenBHadrons_Alt();

      vector<GenBHadron> bcHadrons_alt2;
      if(_isData == false)
	bcHadrons_alt2 = FindGenBHadrons_Alt2();

      // Find Jets
      JetRecord CaloJets = FindCaloJets(bcHadrons,bcHadrons_alt,bcHadrons_alt2,DefaultCaloThreshold);
      //JetRecord MatchedCaloJets = FindMatchedCaloJets(BHadrons, DefaultCaloThreshold);
      //JetRecord UncorrectedCaloJets = FindUncorrectedCaloJets(DefaultUncorrectedCaloThreshold);
      JetRecord PFJets = FindPFJets(TheMuon,bcHadrons, DefaultPFThreshold);
      //JetRecord TrackJets = FindTrackJets(TheMuon, GoodPV, DefaultTrackThreshold);

      /*
      cout << "returned: " << endl;
      for(int j =0 ;j<CaloJets.Count;++j)
	cout << CaloJets.B_DR[j] << " " << CaloJets.B_DR_alt[j] << " " << CaloJets.B_DR_alt2[j] << endl;
      */
      // B-tagging
      int AK5JetCount = nAK5Jet;
      if(AK5JetCount > 200)
	AK5JetCount = 200;

      vector<float> AllBTags(combinedSecondaryVertexBJetTagsAK5Jet,
			     combinedSecondaryVertexBJetTagsAK5Jet + nAK5Jet);
      
      for(int i = 0; i < nAK5Jet; i++)   // zero-out jets out of acceptance.
	{
	  if(sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]) < 30)
            AllBTags[i] = 0;
	  if(etaAK5Jet[i] < -2.4 || etaAK5Jet[i] > 2.4)
            AllBTags[i] = 0;
	}
      
      sort(AllBTags.begin(), AllBTags.end(), greater<float>());

      double HighestBTag = 0;
      double SecondHighestBTag = 0;

      if(AllBTags.size() > 0)
	HighestBTag = AllBTags[0];
      if(AllBTags.size() > 1)
	SecondHighestBTag = AllBTags[1];
      
      // MET-related
      TLorentzVector MuonVector(TheMuon.Px, TheMuon.Py, TheMuon.Pz,
				sqrt(TheMuon.Px * TheMuon.Px + TheMuon.Py * TheMuon.Py + TheMuon.Pz * TheMuon.Pz));
      TLorentzVector METVector(pxMet[0], pyMet[0], pzMet[0], energyMet[0]);
      TLorentzVector PFMETVector(pxPFMet[0], pyPFMet[0], pzPFMet[0], energyPFMet[0]);
      TLorentzVector TCMETVector(pxTCMet[0], pyTCMet[0], pzTCMet[0], energyTCMet[0]);

      // MT's
      TLorentzVector MuonTransverseVector(TheMuon.Px, TheMuon.Py, 0, sqrt(TheMuon.Px * TheMuon.Px + TheMuon.Py * TheMuon.Py));

      TLorentzVector CaloNeutrinoVector = METVector - MuonVector;
      TLorentzVector CaloNeutrinoTransverseVector(CaloNeutrinoVector.Px(), CaloNeutrinoVector.Py(), 0,
						  sqrt(CaloNeutrinoVector.Px() * CaloNeutrinoVector.Px() + CaloNeutrinoVector.Py() * CaloNeutrinoVector.Py()));
      double CaloMT = (CaloNeutrinoTransverseVector + MuonTransverseVector).M();

      TLorentzVector PFNeutrinoVector = PFMETVector;
      TLorentzVector PFNeutrinoTransverseVector(PFNeutrinoVector.Px(), PFNeutrinoVector.Py(), 0,
						sqrt(PFNeutrinoVector.Px() * PFNeutrinoVector.Px() + PFNeutrinoVector.Py() * PFNeutrinoVector.Py()));
      double PFMT = (PFNeutrinoTransverseVector + MuonTransverseVector).M();

      TLorentzVector TCNeutrinoVector = (TCMETVector - MuonVector);
      TLorentzVector TCNeutrinoTransverseVector(TCNeutrinoVector.Px(), TCNeutrinoVector.Py(), 0,
						sqrt(TCNeutrinoVector.Px() * TCNeutrinoVector.Px() + TCNeutrinoVector.Py() * TCNeutrinoVector.Py()));
      double TCMT = (TCNeutrinoTransverseVector + MuonTransverseVector).M();

      // Chris' variables
      double ChrisKinematicTopVariable = 0;
      double MR = 0;
      double MRPrime = 0;
      double R = 0;
      double RPrime = 0;

      // Is there a b-hadron in gen particle list?
      bool HasB = ExistBInGenParticle(bcHadrons);

      // Fill tree
      Messenger.RunNumber = runNumber;
      Messenger.EventNumber = eventNumber;
      Messenger.BunchCrossing = 0;
      Messenger.LumiSection = lumiBlock;
      Messenger.Orbit = orbitNumber;

      Messenger.TheMuon = TheMuon;
      
      // cout << TheMuon.EcalIsolation << endl;
      // cout << TheMuon.HcalIsolation << endl;
      // cout << TheMuon.TrackIsolation << endl;

      Messenger.EcalRelativeIsolation = TheMuon.EcalIsolation / TheMuon.PT;
      if(PFMT > 1e-5)
	Messenger.EcalIsolationOverMT = TheMuon.EcalIsolation / PFMT;
      else
	Messenger.EcalIsolationOverMT = -1;
      
      Messenger.HcalRelativeIsolation = TheMuon.HcalIsolation / TheMuon.PT;
      if(PFMT > 1e-5)
	Messenger.HcalIsolationOverMT = TheMuon.HcalIsolation / PFMT;
      else
	Messenger.HcalIsolationOverMT = -1;
      
      Messenger.TrackRelativeIsolation = TheMuon.TrackIsolation / TheMuon.PT;
      if(PFMT > 1e-5)
	Messenger.TrackIsolationOverMT = TheMuon.TrackIsolation / PFMT;
      else
	Messenger.TrackIsolationOverMT = -1;
      /*
      cout << "before insertion" << endl;
      for(int j =0 ;j<CaloJets.Count;++j)
	cout << CaloJets.B_DR[j] << " " << CaloJets.B_DR_alt[j] << " " << CaloJets.B_DR_alt2[j] << endl;
      */
      Messenger.CaloJets = CaloJets;
      /*
      cout << "after insertion" << endl;
      for(int j = 0;j<Messenger.CaloJets.Count;j++)
	{
	  cout << Messenger.CaloJets.PT[j] << " " << Messenger.CaloJets.Eta[j] << " " << Messenger.CaloJets.B_DR[j] << " " << Messenger.CaloJets.B_DR_alt[j] << " " << Messenger.CaloJets.B_DR_alt2[j] << endl;
	}
      */
      JetRecord UncorrectedCaloJets;
      Messenger.UncorrectedCaloJets = UncorrectedCaloJets;
      JetRecord MatchedCaloJets;
      Messenger.MatchedCaloJets = MatchedCaloJets;
      Messenger.PFJets = PFJets;
      JetRecord TrackJets;
      Messenger.TrackJets = TrackJets;

      Messenger.HighestBTag = HighestBTag;
      Messenger.SecondHighestBTag = SecondHighestBTag;

      Messenger.ChrisKinematicTopVariable = ChrisKinematicTopVariable;
      Messenger.MR = MR;
      Messenger.MRPrime = MRPrime;
      Messenger.R = R;
      Messenger.RPrime = RPrime;

      Messenger.CaloMET = (METVector + MuonVector).Pt();
      Messenger.PFMET = PFMETVector.Pt();
      Messenger.TCMET = TCMETVector.Pt();
   
      Messenger.CaloMT = CaloMT;
      Messenger.PFMT = PFMT;
      Messenger.TCMT = TCMT;

      Messenger.Weight = 1;
      Messenger.HasB = HasB;
      /*
      for(int j = 0;j<Messenger.CaloJets.Count;j++)
	{
	  cout << Messenger.CaloJets.PT[j] << " " << Messenger.CaloJets.Eta[j] << " " << Messenger.CaloJets.B_DR[j] << " " << Messenger.CaloJets.B_DR_alt[j] << " " << Messenger.CaloJets.B_DR_alt2[j] << endl;
	}
      */
      outTree->Fill();
    }

  file->cd();
  outTree->Write();

  // Write histograms
  file->cd();
  WriteQMHistograms(QualityMonitoring1DHistograms);
  WriteQMHistograms(QualityMonitoring2DHistograms);
  HEventCount->Write();

  // Delete Histograms
  DeleteQMHistograms(QualityMonitoring1DHistograms);
  DeleteQMHistograms(QualityMonitoring2DHistograms);
  delete HEventCount;

  file->Close();
}
//---------------------------------------------------------------------------
int CreateWJetDataset::BestPV()
{
  // find the highestpT PV
  double maxpT = -9999.;
  for(int i = 0; i < nPV; i++)
    {
      if(SumPtPV[i] > maxpT)
	{
	  iPV = i;
	  maxpT = SumPtPV[i];
	}
    }

  return iPV;
}
//---------------------------------------------------------------------------
double CreateWJetDataset::pTMuon(int i)
{
  return sqrt(pxMuon[i] * pxMuon[i] + pyMuon[i] * pyMuon[i]);
}
//---------------------------------------------------------------------------
double CreateWJetDataset::SumPt(int iMu, int iZ)
{
  double eta0 = etaMuon[iMu];
  double phi0 = phiMuon[iMu];
  double sumPt_tmp = 0;
  for(int i = 0; i < nTrack; i++)
    {
      if(i == trackIndexMuon[iMu])
	continue; // take out the muon

      if(trackValidHitsTrack[i] < 5)
	continue;                                     // minimum number of hits  XXX

      if(fabs(transvImpactParTrack[i] / transvImpactParErrorTrack[i]) > 5.)
	continue;    // track incompatible with the vertex on (x,y)

      if(fabs(PVzPV[BestPV()] - trackVzTrack[i]) > 0.1)
	continue;              // track incompatible with the vertex on z

      TVector3 v(pxTrack[i], pyTrack[i], pzTrack[i]);
      if(sqrt(pow(v.Eta() - eta0, 2.) + pow(v.Phi() - phi0, 2.)) > 0.5)
	continue; // track outside the cone

      if(v.Pt() < 0.500)
	continue;     // minimum pT             

      if(v.Pt() > 500.)
	continue;     // maximum pT             

      sumPt_tmp += v.Pt();
    }

  return sumPt_tmp; 
}
//---------------------------------------------------------------------------
double CreateWJetDataset::DeltaPhi_PiHalf(double phi1, double phi2)
{
  double dp = fabs(DeltaPhi(phi1, phi2));
  if(dp > asin(1.)) 
    dp = asin(1.) * 2. - dp;
  return dp;
}
//---------------------------------------------------------------------------
vector<MuonCandidate> CreateWJetDataset::MakeMuonCandidates()
{
  // Tight muon:
  //    - Global
  //    - PromptTight
  //    - Tracker
  //    - PT 20
  //    - tracker hits > 10
  //    - pixel hits > 0
  //    - transverse impact parameter < 0.2

  vector<MuonCandidate> Muons;
  Muons.reserve(nMuon);

  Utils anaUtils;

  for(int i = 0; i < nMuon; i++)
    {
      int iTrack = trackIndexMuon[i];
      int iGlobalMuonTrack = combinedTrackIndexMuon[i];
      int iSTAMuonTrack = standAloneTrackIndexMuon[i];

      MuonCandidate Candidate;

      Candidate.IsGlobal = (int)anaUtils.muonIdVal(muonIdMuon[i], bits::AllGlobalMuons);
      Candidate.IsPromptTight = (int)anaUtils.muonIdVal(muonIdMuon[i], bits::GlobalMuonPromptTight);
      Candidate.IsTracker = (int)anaUtils.muonIdVal(muonIdMuon[i], bits::AllTrackerMuons);

      int PixelHit = numberOfValidPixelBarrelHitsTrack[iTrack] + numberOfValidPixelEndcapHitsTrack[iTrack];
      int StripHit = numberOfValidStripTIBHitsTrack[iTrack] + numberOfValidStripTOBHitsTrack[iTrack]
	+ numberOfValidStripTIDHitsTrack[iTrack] + numberOfValidStripTECHitsTrack[iTrack];
      int Chi2 = trackNormalizedChi2GlobalMuonTrack[iGlobalMuonTrack];
      int ValidMuonHit = trackValidHitsSTAMuonTrack[iSTAMuonTrack];
      int MuonStations = 10000;

      Candidate.Charge = chargeMuon[i];

      Candidate.Dxy = eleDxyPV(PVxPV[0], PVyPV[0], PVzPV[0],
			       vertexXMuon[i], vertexYMuon[i], vertexZMuon[i],
			       pxMuon[i], pyMuon[i], pzMuon[i]);
      Candidate.Dz = vertexZMuon[i] - PVzPV[0];
      Candidate.EcalIsolation = emEt03Muon[i];
      Candidate.HcalIsolation = hadEt03Muon[i];
      Candidate.TrackIsolation = sumPt03Muon[i];

      Candidate.Eta = etaMuon[i];
      Candidate.Phi = phiMuon[i];
      Candidate.PT = pTMuon(i);
      Candidate.Px = pxMuon[i];
      Candidate.Py = pyMuon[i];
      Candidate.Pz = pzMuon[i];

      Candidate.PassMuonID = true;
      if(Candidate.IsGlobal == false)
	Candidate.PassMuonID = false;
      if(Candidate.IsTracker == false)
	Candidate.PassMuonID = false;
      if(Candidate.IsPromptTight == false)
	Candidate.PassMuonID = false;
      if(PixelHit < 1)
	Candidate.PassMuonID = false;
      if(StripHit + PixelHit < 10)
	Candidate.PassMuonID = false;
      if(Chi2 > 10)
	Candidate.PassMuonID = false;
      if(ValidMuonHit < 1)
	Candidate.PassMuonID = false;
      
      Candidate.PassFirstLeg = true;
      if(Candidate.PT < 20 || abs(Candidate.Eta) > 2.1 || Candidate.PassMuonID == false)
	Candidate.PassFirstLeg = false;

      Candidate.PassSecondLeg = true;
      if(Candidate.PT < 10 || StripHit + PixelHit < 10)
	Candidate.PassSecondLeg = false;

      Candidate.MuonIndex = i;

      Muons.push_back(Candidate);
    }
   
  return Muons;
}
//---------------------------------------------------------------------------
int CreateWJetDataset::NumberPassFirstLeg(vector<MuonCandidate> &Muons)
{
  int count = 0;

  for(unsigned int i = 0; i < Muons.size(); i++)
    {
      if(Muons[i].PassFirstLeg == true)
	count = count + 1;
    }

  return count;
}
//---------------------------------------------------------------------------
int CreateWJetDataset::NumberPassSecondLeg(vector<MuonCandidate> &Muons)
{
  int count = 0;

  for(unsigned int i = 0; i < Muons.size(); i++)
    {
      if(Muons[i].PassSecondLeg == true)
	count = count + 1;
    }

  return count;
}
//---------------------------------------------------------------------------
JetRecord CreateWJetDataset::FindCaloJets(vector<GenBHadron> & bcHadrons,vector<GenBHadron> & bcHadrons_alt,
					  vector<GenBHadron> & bcHadrons_alt2,double Threshold, double OverallScale, 
					  double EtaSlope)
{
  JetRecord CaloJets;

  multimap<double, int, greater<double> > JetPTMap;
  vector<double> EffectiveJetPT;

  double TotalP[3] = {0, 0, 0};
  int BTaggedCount = 0;
  int LooseBTaggedCount = 0;

  for(int i = 0; i < nAK5Jet; i++)
    {
      double JetPT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);

      if(JetPT * (1 + EtaSlope * fabs(etaAK5Jet[i])) * OverallScale < Threshold)   // PT too small!  Ignore.
	continue;
      if(etaAK5Jet[i] > 2.4 || etaAK5Jet[i] < -2.4)   // out of range!  Ignore.
	continue;

      EffectiveJetPT.push_back(JetPT);

      JetPTMap.insert(pair<double, int>(JetPT, i));

      if(combinedSecondaryVertexBJetTagsAK5Jet[i] > 0.7)
	BTaggedCount = BTaggedCount + 1;
      if(combinedSecondaryVertexBJetTagsAK5Jet[i] > 0.4)
	LooseBTaggedCount = LooseBTaggedCount + 1;

      TotalP[0] = TotalP[0] + pxAK5Jet[i];
      TotalP[1] = TotalP[1] + pyAK5Jet[i];
      TotalP[2] = TotalP[2] + pzAK5Jet[i];
    }

  CaloJets.TotalPT = sqrt(TotalP[0] * TotalP[0] + TotalP[1] * TotalP[1]);

  CaloJets.Count = JetPTMap.size();
  CaloJets.BTaggedCount = BTaggedCount;
  CaloJets.LooseBTaggedCount = LooseBTaggedCount;

  int count = 0;
  for(multimap<double, int, greater<double> >::iterator iter = JetPTMap.begin(); iter != JetPTMap.end(); iter++)
    {
      int index = iter->second;

      CaloJets.PT[count] = sqrt(pxAK5Jet[index] * pxAK5Jet[index] + pyAK5Jet[index] * pyAK5Jet[index]);
      CaloJets.Eta[count] = etaAK5Jet[index];
      CaloJets.Phi[count] = phiAK5Jet[index];
      CaloJets.Energy[count] = energyAK5Jet[index];
      CaloJets.CSV[count] = combinedSecondaryVertexBJetTagsAK5Jet[index];
      CaloJets.CSVMVA[count] = combinedSecondaryVertexMVABJetTagsAK5Jet[index];
      CaloJets.JBP[count] = jetBProbabilityBJetTagsAK5Jet[index];
      CaloJets.JP[count] = jetProbabilityBJetTagsAK5Jet[index];
      CaloJets.TCHE[count] = trackCountingHighEffBJetTagsAK5Jet[index];
      CaloJets.TCHP[count] = trackCountingHighPurBJetTagsAK5Jet[index];
      CaloJets.SSV[count] = simpleSecondaryVertexHighEffBJetTagsAK5Jet[index];



      {
	double minCDR = 9.99;
	double minBDR = 9.99;
	TVector3 JET,HAD;
	JET.SetPtEtaPhi(1,CaloJets.Eta[count],CaloJets.Phi[count]);
	//cout << "bcHadrons " << bcHadrons.size() << endl;
	for(int i =0;i<bcHadrons.size();i++)
	  {
	    int flav = fabs(bcHadrons[i].PDGID);
	    HAD.SetPtEtaPhi(1,bcHadrons[i].Eta,bcHadrons[i].Phi);
	    double dR = HAD.DeltaR(JET);
	    if(flav == 4) //for c hadrons
	      {
		if(minCDR > dR)
		  minCDR = dR;
	      }
	    else if(flav ==5)
	      {
		if(minBDR > dR)
		  minBDR = dR;
	      }
	  }
	
	CaloJets.B_DR[count] = minBDR;
	CaloJets.C_DR[count] = minCDR;
      }

      {
	double minCDR = 9.99;
	double minBDR = 9.99;
	TVector3 JET,HAD;
	JET.SetPtEtaPhi(1,CaloJets.Eta[count],CaloJets.Phi[count]);
	//cout << "bcHadrons_alt " << bcHadrons_alt.size() << endl;
	//cout << "JET: " << JET.Pt() << " " << JET.Eta() << " " << JET.Phi() << endl;
	for(int i =0;i<bcHadrons_alt.size();i++)
	  {
	    int flav = fabs(bcHadrons_alt[i].PDGID);
	    HAD.SetPtEtaPhi(1,bcHadrons_alt[i].Eta,bcHadrons_alt[i].Phi);
	    //cout << "BHAD: " << flav << " " << HAD.Pt() << " " << HAD.Eta() << " " << HAD.Phi() << endl;
	    double dR = HAD.DeltaR(JET);
	    if(flav == 4) //for c hadrons
	      {
		if(minCDR > dR)
		  minCDR = dR;
	      }
	    else if(flav ==5)
	      {
		if(minBDR > dR)
		  minBDR = dR;
	      }
	    //cout << dR << " " << minBDR << " " << minCDR << " " << endl;
	  }
	
	CaloJets.B_DR_alt[count] = minBDR;
	CaloJets.C_DR_alt[count] = minCDR;
	//cout  << CaloJets.B_DR_alt[count] << " " << CaloJets.C_DR_alt[count] << endl;
      }

      {
	double minCDR = 9.99;
	double minBDR = 9.99;
	TVector3 JET,HAD;
	JET.SetPtEtaPhi(1,CaloJets.Eta[count],CaloJets.Phi[count]);
	//cout << "bcHadrons_alt2 " << bcHadrons_alt2.size() << endl;
	for(int i =0;i<bcHadrons_alt2.size();i++)
	  {
	    int flav = fabs(bcHadrons_alt2[i].PDGID);
	    HAD.SetPtEtaPhi(1,bcHadrons_alt2[i].Eta,bcHadrons_alt2[i].Phi);
	    double dR = HAD.DeltaR(JET);
	    if(flav == 4) //for c hadrons
	      {
		if(minCDR > dR)
		  minCDR = dR;
	      }
	    else if(flav ==5)
	      {
		if(minBDR > dR)
		  minBDR = dR;
	      }
	  }
	
	CaloJets.B_DR_alt2[count] = minBDR;
	CaloJets.C_DR_alt2[count] = minCDR;
      }

      //cout << "in loop" << endl;
      //cout << CaloJets.B_DR[count] << " " << CaloJets.B_DR_alt[count] << " " << CaloJets.B_DR_alt2[count] << endl;
      
      count = count + 1;
      
      //change 10 to 20
      if(count >= 20)
	break;
    }
   
  for(unsigned int i = 0; i < EffectiveJetPT.size(); i++)
    {
      int index = (int)EffectiveJetPT[i];
      if(index > 100)
	index = 100;
      if(index < 0)   // huh?
	continue;
      CaloJets.CountSurvey[index] = CaloJets.CountSurvey[index] + 1;
    }
  for(int i = 99; i >= 0; i--)
    CaloJets.CountSurvey[i] = CaloJets.CountSurvey[i] + CaloJets.CountSurvey[i+1];

  return CaloJets;
}

//---------------------------------------------------------------------------
/*
JetRecord CreateWJetDataset::FindUncorrectedCaloJets(double Threshold, double OverallScale, double EtaSlope)
{
  JetRecord UncorrectedCaloJets;

  multimap<double, int, greater<double> > JetPTMap;
  vector<double> EffectiveJetPT;

  double TotalP[3];
  int BTaggedCount = 0;
  int LooseBTaggedCount = 0;

  for(int i = 0; i < nAK5Jet; i++)
    {
      double JetPT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]) / energyAK5Jet[i] * uncorrEnergyAK5Jet[i];
      
      EffectiveJetPT.push_back(JetPT * (1 + EtaSlope * fabs(etaAK5Jet[i])) * OverallScale);

      if(JetPT * (1 + EtaSlope * fabs(etaAK5Jet[i])) * OverallScale < Threshold)   // PT too small!  Ignore.
	continue;
      if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)   // out of range!  Ignore.
	continue;

      JetPTMap.insert(pair<double, int>(JetPT, i));

      if(combinedSecondaryVertexBJetTagsAK5Jet[i] > 0.7)
	BTaggedCount = BTaggedCount + 1;
      if(combinedSecondaryVertexBJetTagsAK5Jet[i] > 0.4)
	LooseBTaggedCount = LooseBTaggedCount + 1;

      TotalP[0] = TotalP[0] + pxAK5Jet[i];
      TotalP[1] = TotalP[1] + pyAK5Jet[i];
      TotalP[2] = TotalP[2] + pzAK5Jet[i];
    }

  UncorrectedCaloJets.TotalPT = sqrt(TotalP[0] * TotalP[0] + TotalP[1] * TotalP[1]);

  UncorrectedCaloJets.Count = JetPTMap.size();
  UncorrectedCaloJets.BTaggedCount = BTaggedCount;
  UncorrectedCaloJets.LooseBTaggedCount = LooseBTaggedCount;

  int count = 0;
  for(multimap<double, int, greater<double> >::iterator iter = JetPTMap.begin(); iter != JetPTMap.end(); iter++)
    {
      int index = iter->second;

      UncorrectedCaloJets.PT[count] = sqrt(pxAK5Jet[index] * pxAK5Jet[index] + pyAK5Jet[index] * pyAK5Jet[index]) / energyAK5Jet[index] * uncorrEnergyAK5Jet[index];
      UncorrectedCaloJets.Eta[count] = etaAK5Jet[index];
      UncorrectedCaloJets.Phi[count] = phiAK5Jet[index];

      count = count + 1;

      if(count >= 10)
	break;
    }

  for(unsigned int i = 0; i < EffectiveJetPT.size(); i++)
    {
      int index = (int)EffectiveJetPT[i];
      if(index > 100)
	index = 100;
      if(index < 0)
	continue;
      UncorrectedCaloJets.CountSurvey[index] = UncorrectedCaloJets.CountSurvey[index] + 1;
    }
  for(int i = 99; i >= 0; i--)
    UncorrectedCaloJets.CountSurvey[i] = UncorrectedCaloJets.CountSurvey[i]
      + UncorrectedCaloJets.CountSurvey[i+1];
   
  return UncorrectedCaloJets;
}
*/
//---------------------------------------------------------------------------
JetRecord CreateWJetDataset::FindPFJets(MuonCandidate &TheMuon,vector<GenBHadron> & bcHadrons, double Threshold, double OverallScale, double EtaSlope)
{
  JetRecord PFJets;

  multimap<double, int, greater<double> > JetPTMap;
  vector<double> EffectiveJetPT;

  double TotalP[3] = {0};
  int BTaggedCount = 0;
  int LooseBTaggedCount = 0;

  for(int j = 0; j < nAK5PFJet; j++)
    {
      // if the jet is muon, don't count
      bool NotMuonJet = true;

      double DEta = etaAK5PFJet[j] - TheMuon.Eta;
      double DPhi = DeltaPhi((double)phiAK5PFJet[j], TheMuon.Phi);
      double Distance = sqrt(DEta * DEta + DPhi * DPhi);

      if(Distance < 0.3)   // note that the isolation cone for muon is 0.3 in VBTF (0.5 in EWK-08-006)
	NotMuonJet = false;

      if(NotMuonJet == false)
	continue;

      double JetPT = sqrt(pxAK5PFJet[j] * pxAK5PFJet[j] + pyAK5PFJet[j] * pyAK5PFJet[j]);

      EffectiveJetPT.push_back(JetPT * (1 + EtaSlope * fabs(etaAK5PFJet[j])) * OverallScale);
      
      if(JetPT * (1 + EtaSlope * fabs(etaAK5PFJet[j])) * OverallScale < Threshold)   // PT too small!  Ignore.
	continue;
      if(etaAK5PFJet[j] > 3 || etaAK5PFJet[j] < -3)   // out of range!  Ignore.
	continue;

      JetPTMap.insert(pair<double, int>(JetPT, j));

      if(combinedSecondaryVertexBJetTagsAK5PFJet[j] > 0.7)
	BTaggedCount = BTaggedCount + 1;
      if(combinedSecondaryVertexBJetTagsAK5PFJet[j] > 0.4)
	LooseBTaggedCount = LooseBTaggedCount + 1;

      TotalP[0] = TotalP[0] + pxAK5PFJet[j];
      TotalP[1] = TotalP[1] + pyAK5PFJet[j];
      TotalP[2] = TotalP[2] + pzAK5PFJet[j];
    }

  PFJets.Count = JetPTMap.size();
  PFJets.BTaggedCount = BTaggedCount;
  PFJets.LooseBTaggedCount = LooseBTaggedCount;

  PFJets.TotalPT = sqrt(TotalP[0] * TotalP[0] + TotalP[1] * TotalP[1]);

  int count = 0;
  for(multimap<double, int, greater<double> >::iterator iter = JetPTMap.begin(); iter != JetPTMap.end(); iter++)
    {
      int index = iter->second;

      PFJets.PT[count] = sqrt(pxAK5PFJet[index] * pxAK5PFJet[index] + pyAK5PFJet[index] * pyAK5PFJet[index]);
      PFJets.Eta[count] = etaAK5PFJet[index];
      PFJets.Phi[count] = phiAK5PFJet[index];
      PFJets.Energy[count] = energyAK5PFJet[index];
      PFJets.CSV[count] = combinedSecondaryVertexBJetTagsAK5PFJet[index];
      PFJets.CSVMVA[count] = combinedSecondaryVertexMVABJetTagsAK5PFJet[index];
      PFJets.JBP[count] = jetBProbabilityBJetTagsAK5PFJet[index];
      PFJets.JP[count] = jetProbabilityBJetTagsAK5PFJet[index];
      PFJets.TCHE[count] = trackCountingHighEffBJetTagsAK5PFJet[index];
      PFJets.TCHP[count] = trackCountingHighPurBJetTagsAK5PFJet[index];
      PFJets.SSV[count] = simpleSecondaryVertexHighEffBJetTagsAK5PFJet[index];
      
      double minCDR = 9.99;
      double minBDR = 9.99;
      TVector3 JET,HAD;
      JET.SetPtEtaPhi(1,PFJets.Eta[count],PFJets.Phi[count]);
      for(int i =0;i<bcHadrons.size();i++)
	{
	  int flav = fabs(bcHadrons[i].PDGID);
	  HAD.SetPtEtaPhi(1,bcHadrons[i].Eta,bcHadrons[i].Phi);
	  double dR = HAD.DeltaR(JET);
	  if(flav == 4) //for c hadrons
	    {
	      if(minCDR > dR)
		minCDR = dR;
	    }
	  else if(flav ==5)
	    {
	      if(minBDR > dR)
		minBDR = dR;
	    }
	}

      PFJets.B_DR[count] = minBDR;
      PFJets.C_DR[count] = minCDR;

      count = count + 1;

      if(count >= 20)
	break;
    }
   
  for(unsigned int i = 0; i < EffectiveJetPT.size(); i++)
    {
      int index = (int)EffectiveJetPT[i];
      if(index > 100)
	index = 100;
      if(index < 0)
	continue;
      PFJets.CountSurvey[index] = PFJets.CountSurvey[index] + 1;
    }
  for(int i = 99; i >= 0; i--)
    PFJets.CountSurvey[i] = PFJets.CountSurvey[i] + PFJets.CountSurvey[i+1];

  return PFJets;
}
//---------------------------------------------------------------------------
/*
JetRecord CreateWJetDataset::FindTrackJets(MuonCandidate &TheMuon, int GoodPV, double Threshold,
					   double OverallScale, double EtaSlope)
{
  JetRecord TrackJetRecord;

  // candidate muon
  int CandidateMuonTrackIndex = -1;
  if(TheMuon.MuonIndex >= 0 && TheMuon.MuonIndex < nTrack)
    CandidateMuonTrackIndex = trackIndexMuon[TheMuon.MuonIndex];

  // Prepare lorentz vectors for the components - remove candidate muons
  vector<TLorentzVector> tracks;
  for(int j = 0; j < nTrack; j++)
    {
      if(CandidateMuonTrackIndex == j)
	continue;

      TVector3 TrackVertex(trackVxTrack[j], trackVyTrack[j], trackVzTrack[j]);
      TVector3 PrimaryVertex(PVxPV[GoodPV], PVyPV[GoodPV], PVzPV[GoodPV]);
      TrackVertex = TrackVertex - PrimaryVertex;

      TVector3 TrackMomentum(pxTrack[j], pyTrack[j], pzTrack[j]);

      // basic track quality
      if(trackValidHitsTrack[j] < 5)
	continue;
      if(trackNormalizedChi2Track[j] > 20)
	continue;
      
      // consistent with primary vertex - z direction 0.1 cm, xy direction 600 mu m
      if(fabs(TrackVertex.Z()) > 0.1)
	continue;
      if(fabs(TrackVertex.Mag()) > 0.1)
	continue;
      if(fabs(transvImpactParTrack[j]) > 0.06)
	continue;

      // there should be no closer vertex
      bool CloserVertexFound = false;
      for(int k = 0; k < nPV; k++)
	{
	  if(k == iPV)
            continue;
      
	  if(fabs(trackVzTrack[j] - PVzPV[GoodPV]) > fabs(trackVzTrack[j] - PVzPV[k]))
            CloserVertexFound = true;
	}
      
      if(CloserVertexFound == true)
	continue;

      // track vertex should be close to muon candidate vertex (To be added.  Will be useful in pileup.)

      // track pt and eta requirement
      if(TrackMomentum.Pt() < 0.5 || TrackMomentum.Pt() > 500)
	continue;
      if(TrackMomentum.Eta() < -2.4 || TrackMomentum.Eta() > 2.4)
	continue;

      TLorentzVector FourMomentum;
      FourMomentum.SetPxPyPzE(pxTrack[j], pyTrack[j], pzTrack[j], TrackMomentum.Mag());

      tracks.push_back(FourMomentum);
    }

  if(tracks.size() == 0)
    {
      TrackJetRecord.Count = 0;
      return TrackJetRecord;
    }

  // do clustering - what is the second parameter?   Ans: pt_min
  vector<Jet> TrackJets = FastJetAlgorithm(tracks, 0.5, 0.0);

  // export number of jets and also the highest 3 jets (in pt)
  multimap<double, int, greater<double> > JetPTMap;
  vector<double> EffectiveJetPT;

  double TotalP[3];

  for(unsigned int j = 0; j < TrackJets.size(); j++)
    {
      EffectiveJetPT.push_back(TrackJets[j].Pt() * (1 + EtaSlope * fabs(TrackJets[j].Eta())) * OverallScale);
      
      if(TrackJets[j].Pt() * (1 + EtaSlope * fabs(TrackJets[j].Eta())) * OverallScale < Threshold)
	continue;
      if(TrackJets[j].Eta() > 2.4 || TrackJets[j].Eta() < -2.4)
	continue;

      JetPTMap.insert(pair<double, int>(TrackJets[j].Pt(), j));

      TotalP[0] = TotalP[0] + TrackJets[j].Px();
      TotalP[1] = TotalP[1] + TrackJets[j].Py();
      TotalP[2] = TotalP[2] + TrackJets[j].Pz();
    }

  TrackJetRecord.TotalPT = sqrt(TotalP[0] * TotalP[0] + TotalP[1] * TotalP[1]);

  TrackJetRecord.Count = JetPTMap.size();
  TrackJetRecord.BTaggedCount = 0;   // I don't have time to re-implement b-tagging here....
  TrackJetRecord.LooseBTaggedCount = 0;   // let's fill 0 for now

  int count = 0;
  for(multimap<double, int>::iterator iter = JetPTMap.begin(); iter != JetPTMap.end(); iter++)
    {
      int index = iter->second;

      TrackJetRecord.PT[count] = TrackJets[index].Pt();
      TrackJetRecord.Eta[count] = TrackJets[index].Eta();
      TrackJetRecord.Phi[count] = TrackJets[index].Phi();

      count = count + 1;

      if(count >= 10)
	break;
    }

  for(unsigned int i = 0; i < EffectiveJetPT.size(); i++)
    {
      int index = (int)EffectiveJetPT[i];
      if(index > 100)
	index = 100;
      if(index < 0)
	continue;
      TrackJetRecord.CountSurvey[index] = TrackJetRecord.CountSurvey[index] + 1;
    }
  for(int i = 99; i >= 0; i--)
    TrackJetRecord.CountSurvey[i] = TrackJetRecord.CountSurvey[i] + TrackJetRecord.CountSurvey[i+1];

  return TrackJetRecord;
}

*/
vector<GenBHadron> CreateWJetDataset::FindGenBHadrons_Alt()
{
   vector<GenBHadron> BHadrons;

   for(int i = 0; i < nMc; i++)
   {
     //cout << "- " << i << "\t" << statusMc[i] << " " << idMc[i] << endl;
     if(! (statusMc[i]==3 && (fabs(idMc[i])==4 || fabs(idMc[i])==5)))   // We want b and c hadrons
       continue;
     //cout << "accepted" << endl;
     // if(statusMc[i] != 2)   // status-2 only
     //    continue;
     
     GenBHadron HadronCandidate;
     
     HadronCandidate.Final = true;
     HadronCandidate.Eta = etaMc[i];
     HadronCandidate.Phi = phiMc[i];
     HadronCandidate.Energy = energyMc[i];
     HadronCandidate.Status = statusMc[i];
     HadronCandidate.PDGID = abs(idMc[i]);
     
     //VecbosIndexToVectorIndex.insert(pair<int, int>(i, BHadrons.size()));
     //VectorIndexToVecbosIndex.push_back(i);
     
     BHadrons.push_back(HadronCandidate);
   }

   /*
   for(int i = 0; i < (int)BHadrons.size(); i++)
     {
       int MotherVecbosIndex = mothMc[VectorIndexToVecbosIndex[i]];
       if(VecbosIndexToVectorIndex.find(MotherVecbosIndex) == VecbosIndexToVectorIndex.end())   // mother not b-hadron
         continue;
       
       BHadrons[VecbosIndexToVectorIndex[MotherVecbosIndex]].Final = false;
     }
   */
   // int BHadronCount = 0;
   // for(int i = 0; i < (int)BHadrons.size(); i++)
   //    if(BHadrons[i].Final == true)
   //       BHadronCount = BHadronCount + 1;
   // cout << "Number of final b-hadrons: " << BHadronCount << "/" << BHadrons.size() << endl;

   return BHadrons;
}


vector<GenBHadron> CreateWJetDataset::FindGenBHadrons_Alt2()
{
  vector<GenBHadron> BHadrons;
  
  //cout << "LIST ALL MC's, STORE DAUGHTERS" << endl;
  map<int,vector<int> > daughters;
  for(int i = 0; i < nMc; i++)
    {
      //cout << i << " " << idMc[i] << " " << mothMc[i];
      daughters[mothMc[i]].push_back(i);
    }
    
  
  //cout << "FILTER MOST FINAL B AND C HADRONS" << endl;
  for(int i = 0; i < nMc; i++)
    {
      //cout << i << " " << idMc[i] << " ";
      //if(daughters.find(i)!=daughters.end())
      //{
      //  cout << "d: ";
      //  for(int d = 0;d<daughters[i].size();++d)
      //    cout << daughters[i][d] << " " << idMc[daughters[i][d]] << " " << mothMc[daughters[i][d]] << " || ";
      //  cout << endl;
      //}
      //else
      //cout << endl;

      int p_u = abs(idMc[i]);
      int p_h = abs(idMc[i]%1000/100);
      int p_t = abs(idMc[i]%10000/1000);
      int taste = -1;
      //if the particle is a b meson or baryon                                                                                                             
      if(p_h==5 || p_t == 5 || p_u ==5)
	taste = 5;
      //if the particle is a c meson or baryon                                                                                                             
      else if(p_h==4 || p_t == 4 || p_u ==4)
	taste = 4;
      else
	continue;


      bool store = false;
      if(taste >0)
	{
	  //store the particle                                                                                                                             
	  store = true;
	  //if it decays to a particle of different 'taste'                                                                                                
	  for(unsigned int d = 0;d<daughters[i].size();++d)
	    {
	      int p_u = abs(idMc[daughters[i][d]]);
	      int d_h = abs(idMc[daughters[i][d]]%1000/100);
	      int d_t = abs(idMc[daughters[i][d]]%10000/1000);
	      if(d_h==taste||d_t==taste || p_u == taste)
		store = false;
	    }
	}

      if(!store)
	continue;

      
      
      GenBHadron HadronCandidate;
      
      HadronCandidate.Final = true;
      HadronCandidate.Eta = etaMc[i];
      HadronCandidate.Phi = phiMc[i];
      HadronCandidate.Energy = energyMc[i];
      HadronCandidate.Status = statusMc[i];
      HadronCandidate.PDGID = taste;
      
      //VecbosIndexToVectorIndex.insert(pair<int, int>(i, BHadrons.size()));
      //VectorIndexToVecbosIndex.push_back(i);
      
      BHadrons.push_back(HadronCandidate);
    }
  
   /*
   for(int i = 0; i < (int)BHadrons.size(); i++)
     {
       int MotherVecbosIndex = mothMc[VectorIndexToVecbosIndex[i]];
       if(VecbosIndexToVectorIndex.find(MotherVecbosIndex) == VecbosIndexToVectorIndex.end())   // mother not b-hadron
         continue;
       
       BHadrons[VecbosIndexToVectorIndex[MotherVecbosIndex]].Final = false;
     }
   */
   // int BHadronCount = 0;
   // for(int i = 0; i < (int)BHadrons.size(); i++)
   //    if(BHadrons[i].Final == true)
   //       BHadronCount = BHadronCount + 1;
   // cout << "Number of final b-hadrons: " << BHadronCount << "/" << BHadrons.size() << endl;

   return BHadrons;
}


//---------------------------------------------------------------------------
/*
bool CreateWJetDataset::IsBHadron(int PDGID)
{
   if(PDGID < 0)
      PDGID = -PDGID;

   if(PDGID == 5)
      return true;
   if((PDGID % 1000) / 100 == 5)
      return true;
   if((PDGID % 10000) / 1000 == 5)
      return true;

   return false;
}
*/
//---------------------------------------------------------------------------
void CreateWJetDataset::MakeBranches(TTree *tree, OutputRecord *record)
{
  // General
  tree->Branch("RunNumber", &record->RunNumber, "RunNumber/D");
  tree->Branch("EventNumber", &record->EventNumber, "EventNumber/D");
  tree->Branch("BunchCrossing", &record->BunchCrossing, "BunchCrossing/D");
  tree->Branch("LumiSection", &record->LumiSection, "LumiSection/D");
  tree->Branch("Orbit", &record->Orbit, "Orbit/D");

  tree->Branch("PassHLT", &record->PassHLT, "PassHLT/O");

  // The muon
  tree->Branch("Muon1PT", &record->TheMuon.PT, "Muon1PT/D");
  tree->Branch("Muon1Eta", &record->TheMuon.Eta, "Muon1Eta/D");
  tree->Branch("Muon1Phi", &record->TheMuon.Phi, "Muon1Phi/D");
  tree->Branch("Muon1PassMuonID", &record->TheMuon.PassMuonID, "Muon1PassMuonID/O");
  tree->Branch("Muon1PassFirstLeg", &record->TheMuon.PassFirstLeg, "Muon1PassFirstLeg/O");
  tree->Branch("Muon1PassSecondLeg", &record->TheMuon.PassSecondLeg, "Muon1PassSecondLeg/O");
  tree->Branch("Muon1Dxy", &record->TheMuon.Dxy, "Muon1Dxy/D");
  tree->Branch("Muon1Dz", &record->TheMuon.Dz, "Muon1Dz/D");
  tree->Branch("Muon1Charge", &record->TheMuon.Charge, "Muon1Charge/I");
  tree->Branch("Muon1EcalRelativeIsolation", &record->EcalRelativeIsolation, "Muon1EcalRelativeIsolation/D");
  tree->Branch("Muon1EcalIsolationOverMT", &record->EcalIsolationOverMT, "Muon1EcalIsolationOverMT/D");
  tree->Branch("Muon1HcalRelativeIsolation", &record->HcalRelativeIsolation, "Muon1HcalRelativeIsolation/D");
  tree->Branch("Muon1HcalIsolationOverMT", &record->HcalIsolationOverMT, "Muon1HcalIsolationOverMT/D");
  tree->Branch("Muon1TrackRelativeIsolation", &record->TrackRelativeIsolation, "Muon1TrackRelativeIsolation/D");
  tree->Branch("Muon1TrackIsolationOverMT", &record->TrackIsolationOverMT, "Muon1TrackIsolationOverMT/D");

  // Jets
  tree->Branch("NCaloJet", &record->CaloJets.Count, "NCaloJet/I");
  tree->Branch("NBTaggedCaloJet", &record->CaloJets.BTaggedCount, "NBTaggedCaloJet/I");
  tree->Branch("NLooseBTaggedCaloJet", &record->CaloJets.LooseBTaggedCount, "NLooseBTaggedCaloJet/I");
  tree->Branch("CaloJet1PT", &record->CaloJets.PT[0], "CaloJet1PT/D");
  tree->Branch("CaloJet1Eta", &record->CaloJets.Eta[0], "CaloJet1Eta/D");
  tree->Branch("CaloJet1Phi", &record->CaloJets.Phi[0], "CaloJet1Phi/D");
  tree->Branch("CaloJet2PT", &record->CaloJets.PT[1], "CaloJet2PT/D");
  tree->Branch("CaloJet2Eta", &record->CaloJets.Eta[1], "CaloJet2Eta/D");
  tree->Branch("CaloJet2Phi", &record->CaloJets.Phi[1], "CaloJet2Phi/D");
  tree->Branch("CaloJet3PT", &record->CaloJets.PT[2], "CaloJet3PT/D");
  tree->Branch("CaloJet3Eta", &record->CaloJets.Eta[2], "CaloJet3Eta/D");
  tree->Branch("CaloJet3Phi", &record->CaloJets.Phi[2], "CaloJet3Phi/D");
  //for(int i = 0; i <= 100; i++)
  //tree->Branch(Form("NCaloJet%d", i), &record->CaloJets.CountSurvey[i], Form("NCaloJet%d/I"));
  /*
  tree->Branch("NMatchedCaloJet", &record->MatchedCaloJets.Count, "NMatchedCaloJet/I");
  tree->Branch("NBTaggedMatchedCaloJet", &record->MatchedCaloJets.BTaggedCount, "NBTaggedMatchedCaloJet/I");
  tree->Branch("NLooseBTaggedMatchedCaloJet", &record->MatchedCaloJets.LooseBTaggedCount, "NLooseBTaggedMatchedCaloJet/I");
  tree->Branch("MatchedCaloJet1PT", &record->MatchedCaloJets.PT[0], "MatchedCaloJet1PT/D");
  tree->Branch("MatchedCaloJet1Eta", &record->MatchedCaloJets.Eta[0], "MatchedCaloJet1Eta/D");
  tree->Branch("MatchedCaloJet1Phi", &record->MatchedCaloJets.Phi[0], "MatchedCaloJet1Phi/D");
  tree->Branch("MatchedCaloJet2PT", &record->MatchedCaloJets.PT[1], "MatchedCaloJet2PT/D");
  tree->Branch("MatchedCaloJet2Eta", &record->MatchedCaloJets.Eta[1], "MatchedCaloJet2Eta/D");
  tree->Branch("MatchedCaloJet2Phi", &record->MatchedCaloJets.Phi[1], "MatchedCaloJet2Phi/D");
  tree->Branch("MatchedCaloJet3PT", &record->MatchedCaloJets.PT[2], "MatchedCaloJet3PT/D");
  tree->Branch("MatchedCaloJet3Eta", &record->MatchedCaloJets.Eta[2], "MatchedCaloJet3Eta/D");
  tree->Branch("MatchedCaloJet3Phi", &record->MatchedCaloJets.Phi[2], "MatchedCaloJet3Phi/D");
  for(int i = 0; i <= 100; i++)
    tree->Branch(Form("NMatchedCaloJet%d", i), &record->MatchedCaloJets.CountSurvey[i], Form("NMatchedCaloJet%d/I"));

  tree->Branch("NUncorrectedCaloJet", &record->UncorrectedCaloJets.Count, "NUncorrectedCaloJet/I");
  tree->Branch("NBTaggedUncorrectedCaloJet", &record->UncorrectedCaloJets.BTaggedCount, "NBTaggedUncorrectedCaloJet/I");
  tree->Branch("NLooseBTaggedUncorrectedCaloJet", &record->UncorrectedCaloJets.LooseBTaggedCount, "NLooseBTaggedUncorrectedCaloJet/I");
  tree->Branch("UncorrectedCaloJet1PT", &record->UncorrectedCaloJets.PT[0], "UncorrectedCaloJet1PT/D");
  tree->Branch("UncorrectedCaloJet1Eta", &record->UncorrectedCaloJets.Eta[0], "UncorrectedCaloJet1Eta/D");
  tree->Branch("UncorrectedCaloJet1Phi", &record->UncorrectedCaloJets.Phi[0], "UncorrectedCaloJet1Phi/D");
  tree->Branch("UncorrectedCaloJet2PT", &record->UncorrectedCaloJets.PT[1], "UncorrectedCaloJet2PT/D");
  tree->Branch("UncorrectedCaloJet2Eta", &record->UncorrectedCaloJets.Eta[1], "UncorrectedCaloJet2Eta/D");
  tree->Branch("UncorrectedCaloJet2Phi", &record->UncorrectedCaloJets.Phi[1], "UncorrectedCaloJet2Phi/D");
  tree->Branch("UncorrectedCaloJet3PT", &record->UncorrectedCaloJets.PT[2], "UncorrectedCaloJet3PT/D");
  tree->Branch("UncorrectedCaloJet3Eta", &record->UncorrectedCaloJets.Eta[2], "UncorrectedCaloJet3Eta/D");
  tree->Branch("UncorrectedCaloJet3Phi", &record->UncorrectedCaloJets.Phi[2], "UncorrectedCaloJet3Phi/D");
  for(int i = 0; i <= 100; i++)
    tree->Branch(Form("NUncorrectedCaloJet%d", i), &record->UncorrectedCaloJets.CountSurvey[i],
		 Form("NUncorrectedCaloJet%d/I"));
  */
  tree->Branch("NPFJet", &record->PFJets.Count, "NPFJet/I");
  tree->Branch("NBTaggedPFJet", &record->PFJets.BTaggedCount, "NBTaggedPFJet/I");
  tree->Branch("NLooseBTaggedPFJet", &record->PFJets.LooseBTaggedCount, "NLooseBTaggedPFJet/I");
  tree->Branch("PFJet1PT", &record->PFJets.PT[0], "PFJet1PT/D");
  tree->Branch("PFJet1Eta", &record->PFJets.Eta[0], "PFJet1Eta/D");
  tree->Branch("PFJet1Phi", &record->PFJets.Phi[0], "PFJet1Phi/D");
  tree->Branch("PFJet2PT", &record->PFJets.PT[1], "PFJet2PT/D");
  tree->Branch("PFJet2Eta", &record->PFJets.Eta[1], "PFJet2Eta/D");
  tree->Branch("PFJet2Phi", &record->PFJets.Phi[1], "PFJet2Phi/D");
  tree->Branch("PFJet3PT", &record->PFJets.PT[2], "PFJet3PT/D");
  tree->Branch("PFJet3Eta", &record->PFJets.Eta[2], "PFJet3Eta/D");
  tree->Branch("PFJet3Phi", &record->PFJets.Phi[2], "PFJet3Phi/D");
  //for(int i = 0; i <= 100; i++)
  //tree->Branch(Form("NPFJet%d", i), &record->PFJets.CountSurvey[i], Form("NPFJet%d/I"));

  /*
  tree->Branch("NTrackJet", &record->TrackJets.Count, "NTrackJet/I");
  tree->Branch("TrackJet1PT", &record->TrackJets.PT[0], "TrackJet1PT/D");
  tree->Branch("TrackJet1Eta", &record->TrackJets.Eta[0], "TrackJet1Eta/D");
  tree->Branch("TrackJet1Phi", &record->TrackJets.Phi[0], "TrackJet1Phi/D");
  tree->Branch("TrackJet2PT", &record->TrackJets.PT[1], "TrackJet2PT/D");
  tree->Branch("TrackJet2Eta", &record->TrackJets.Eta[1], "TrackJet2Eta/D");
  tree->Branch("TrackJet2Phi", &record->TrackJets.Phi[1], "TrackJet2Phi/D");
  tree->Branch("TrackJet3PT", &record->TrackJets.PT[2], "TrackJet3PT/D");
  tree->Branch("TrackJet3Eta", &record->TrackJets.Eta[2], "TrackJet3Eta/D");
  tree->Branch("TrackJet3Phi", &record->TrackJets.Phi[2], "TrackJet3Phi/D");
  for(int i = 0; i <= 100; i++)
    tree->Branch(Form("NTrackJet%d", i), &record->TrackJets.CountSurvey[i], Form("NTrackJet%d/I"));
  */
  // B-Tagging variables
  tree->Branch("HighestBTag", &record->HighestBTag, "HighestBTag/D");
  tree->Branch("SecondHighestBTag", &record->SecondHighestBTag, "SecondHighestBTag/D");

  // Chris variables
  tree->Branch("ChrisKinematicTopVariable", &record->ChrisKinematicTopVariable, "ChrisKinematicTopVariable/D");
  tree->Branch("MR", &record->MR, "MR/D");
  tree->Branch("MRPrime", &record->MRPrime, "MRPrime/D");
  tree->Branch("R", &record->R, "R/D");
  tree->Branch("RPrime", &record->RPrime, "RPrime/D");

  // MET, MT
  tree->Branch("CaloMET", &record->CaloMET, "CaloMET/D");
  tree->Branch("PFMET", &record->PFMET, "PFMET/D");
  tree->Branch("TCMET", &record->TCMET, "TCMET/D");

  tree->Branch("CaloMT", &record->CaloMT, "CaloMT/D");
  tree->Branch("PFMT", &record->PFMT, "PFMT/D");
  tree->Branch("TCMT", &record->TCMT, "TCMT/D");

  // MC-related
  tree->Branch("HasB", &record->HasB, "HasB/O");

  // MC Weight
  tree->Branch("Weight", &record->Weight, "Weight/D");

  //arrays
  tree->Branch("CaloJet_PT",record->CaloJets.PT,"CaloJet_PT[NCaloJet]/D");
  tree->Branch("CaloJet_Eta",record->CaloJets.Eta,"CaloJet_Eta[NCaloJet]/D");
  tree->Branch("CaloJet_Phi",record->CaloJets.Phi,"CaloJet_Phi[NCaloJet]/D");
  tree->Branch("CaloJet_Energy",record->CaloJets.Energy,"CaloJet_Energy[NCaloJet]/D");
  tree->Branch("CaloJet_CSV",record->CaloJets.CSV,"CaloJet_CSV[NCaloJet]/D");
  tree->Branch("CaloJet_CSVMVA",record->CaloJets.CSVMVA,"CaloJet_CSVMVA[NCaloJet]/D");
  tree->Branch("CaloJet_JBP",record->CaloJets.JBP,"CaloJet_JBP[NCaloJet]/D");
  tree->Branch("CaloJet_JP",record->CaloJets.JP,"CaloJet_JP[NCaloJet]/D");
  tree->Branch("CaloJet_TCHE",record->CaloJets.TCHE,"CaloJet_TCHE[NCaloJet]/D");
  tree->Branch("CaloJet_TCHP",record->CaloJets.TCHP,"CaloJet_TCHP[NCaloJet]/D");
  tree->Branch("CaloJet_SSV",record->CaloJets.SSV,"CaloJet_SSV[NCaloJet]/D");
  tree->Branch("CaloJet_B_DR",record->CaloJets.B_DR,"CaloJet_B_DR[NCaloJet]/D");
  tree->Branch("CaloJet_C_DR",record->CaloJets.C_DR,"CaloJet_C_DR[NCaloJet]/D");
  tree->Branch("CaloJet_B_DR_alt",record->CaloJets.B_DR_alt,"CaloJet_B_DR_alt[NCaloJet]/D");
  tree->Branch("CaloJet_C_DR_alt",record->CaloJets.C_DR_alt,"CaloJet_C_DR_alt[NCaloJet]/D");
  tree->Branch("CaloJet_B_DR_alt2",record->CaloJets.B_DR_alt2,"CaloJet_B_DR_alt2[NCaloJet]/D");
  tree->Branch("CaloJet_C_DR_alt2",record->CaloJets.C_DR_alt2,"CaloJet_C_DR_alt2[NCaloJet]/D");


  tree->Branch("PFJet_PT",record->PFJets.PT,"PFJet_PT[NPFJet]/D");
  tree->Branch("PFJet_Eta",record->PFJets.Eta,"PFJet_Eta[NPFJet]/D");
  tree->Branch("PFJet_Phi",record->PFJets.Phi,"PFJet_Phi[NPFJet]/D");
  tree->Branch("PFJet_Energy",record->PFJets.Energy,"PFJet_Energy[NPFJet]/D");
  tree->Branch("PFJet_CSV",record->PFJets.CSV,"PFJet_CSV[NPFJet]/D");
  tree->Branch("PFJet_CSVMVA",record->PFJets.CSVMVA,"PFJet_CSVMVA[NPFJet]/D");
  tree->Branch("PFJet_JBP",record->PFJets.JBP,"PFJet_JBP[NPFJet]/D");
  tree->Branch("PFJet_JP",record->PFJets.JP,"PFJet_JP[NPFJet]/D");
  tree->Branch("PFJet_TCHE",record->PFJets.TCHE,"PFJet_TCHE[NPFJet]/D");
  tree->Branch("PFJet_TCHP",record->PFJets.TCHP,"PFJet_TCHP[NPFJet]/D");
  tree->Branch("PFJet_SSV",record->PFJets.SSV,"PFJet_SSV[NPFJet]/D");
  tree->Branch("PFJet_B_DR",record->PFJets.B_DR,"PFJet_B_DR[NPFJet]/D");
  tree->Branch("PFJet_C_DR",record->PFJets.C_DR,"PFJet_C_DR[NPFJet]/D");
}



//---------------------------------------------------------------------------
map<string, TH1D *> CreateWJetDataset::GenerateQM1DHistograms()
{
  map<string, TH1D *> Histograms;

  return Histograms;
}
//---------------------------------------------------------------------------
map<string, TH2D *> CreateWJetDataset::GenerateQM2DHistograms()
{
  map<string, TH2D *> Histograms;

  return Histograms;
}
//---------------------------------------------------------------------------
void CreateWJetDataset::WriteQMHistograms(map<string, TH1D *> Histograms)
{
  for(map<string, TH1D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
    {
      if(iter->second != NULL)
	iter->second->Write();
    }
}
//---------------------------------------------------------------------------
void CreateWJetDataset::WriteQMHistograms(map<string, TH2D *> Histograms)
{
  for(map<string, TH2D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
    {
      if(iter->second != NULL)
	iter->second->Write();
    }
}
//---------------------------------------------------------------------------
void CreateWJetDataset::DeleteQMHistograms(map<string, TH1D *> Histograms)
{
  for(map<string, TH1D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
    {
      if(iter->second == NULL)
	continue;

      delete iter->second;
      iter->second = NULL;
    }
}
//---------------------------------------------------------------------------
void CreateWJetDataset::DeleteQMHistograms(map<string, TH2D *> Histograms)
{
  for(map<string, TH2D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
    {
      if(iter->second == NULL)
	continue;

      delete iter->second;
      iter->second = NULL;
    }
}
//---------------------------------------------------------------------------
bool CreateWJetDataset::ExistBInGenParticle(vector<GenBHadron> & bcHadrons)
{
  if(_isData == true)
    return false;

  for(int i = 0; i < bcHadrons.size();i++)
    {
      if(fabs(bcHadrons[i].PDGID)==5)
	return true;
    }
  return false;
}
//---------------------------------------------------------------------------













