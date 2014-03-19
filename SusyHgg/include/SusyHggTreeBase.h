//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 14 17:46:02 2014 by ROOT version 5.34/07
// from TTree SusyHggTree/
// found on file: /home/amott/raid4/SMS-TChiHH_2b2g_2J_mChargino-130to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola__Summer12-START53_V19_FSIM-v1_SUSY.root
//////////////////////////////////////////////////////////

#ifndef SusyHggTreeBase_h
#define SusyHggTreeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SusyHggTreeBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Float_t         mgg;
   Float_t         ptgg;
   Float_t         etagg;
   Float_t         phigg;
   Float_t         pho1_pt;
   Float_t         pho1_eta;
   Float_t         pho1_phi;
   Float_t         pho1_r9;
   Float_t         pho1_sigEoE;
   Float_t         pho1_energyGen;
   Char_t          pho1_genMatch;
   Float_t         pho1_sieie;
   Float_t         pho1_HE;
   Float_t         pho1_charged;
   Float_t         pho1_neutral;
   Float_t         pho1_photon;
   Char_t          pho1_eleveto;
   Char_t          pho1_pass_id;
   Char_t          pho1_pass_iso;
   Float_t         pho2_pt;
   Float_t         pho2_eta;
   Float_t         pho2_phi;
   Float_t         pho2_r9;
   Float_t         pho2_sigEoE;
   Float_t         pho2_energyGen;
   Char_t          pho2_genMatch;
   Float_t         pho2_sieie;
   Float_t         pho2_HE;
   Float_t         pho2_charged;
   Float_t         pho2_neutral;
   Float_t         pho2_photon;
   Char_t          pho2_eleveto;
   Char_t          pho2_pass_id;
   Char_t          pho2_pass_iso;
   Float_t         ele1_pt;
   Float_t         ele1_eta;
   Float_t         ele1_phi;
   Float_t         mu1_pt;
   Float_t         mu1_eta;
   Float_t         mu1_phi;
   Float_t         highest_csv;
   Float_t         highest_csv_pt;
   Float_t         highest_csv_eta;
   Float_t         highest_csv_phi;
   Float_t         highest_csv_up;
   Float_t         highest_csv_pt_up;
   Float_t         highest_csv_eta_up;
   Float_t         highest_csv_phi_up;
   Float_t         highest_csv_down;
   Float_t         highest_csv_pt_down;
   Float_t         highest_csv_eta_down;
   Float_t         highest_csv_phi_down;
   Float_t         second_csv;
   Float_t         second_csv_pt;
   Float_t         second_csv_eta;
   Float_t         second_csv_phi;
   Float_t         second_csv_up;
   Float_t         second_csv_pt_up;
   Float_t         second_csv_eta_up;
   Float_t         second_csv_phi_up;
   Float_t         second_csv_down;
   Float_t         second_csv_pt_down;
   Float_t         second_csv_eta_down;
   Float_t         second_csv_phi_down;
   Float_t         mbb;
   Float_t         mbb_up;
   Float_t         mbb_down;
   Float_t         mbb_NearH;
   Float_t         mbb_NearH_up;
   Float_t         mbb_NearH_down;
   Float_t         mbb_NearZ;
   Float_t         mbb_NearZ_up;
   Float_t         mbb_NearZ_down;
   Float_t         hem1_pt;
   Float_t         hem1_eta;
   Float_t         hem1_phi;
   Float_t         hem1_M;
   Float_t         hem1_pt_up;
   Float_t         hem1_eta_up;
   Float_t         hem1_phi_up;
   Float_t         hem1_M_up;
   Float_t         hem1_pt_down;
   Float_t         hem1_eta_down;
   Float_t         hem1_phi_down;
   Float_t         hem1_M_down;
   Float_t         hem2_pt;
   Float_t         hem2_eta;
   Float_t         hem2_phi;
   Float_t         hem2_M;
   Float_t         hem2_pt_up;
   Float_t         hem2_eta_up;
   Float_t         hem2_phi_up;
   Float_t         hem2_M_up;
   Float_t         hem2_pt_down;
   Float_t         hem2_eta_down;
   Float_t         hem2_phi_down;
   Float_t         hem2_M_down;
   Float_t         MET;
   Float_t         METPhi;
   Int_t           Njets;
   Float_t         MR;
   Float_t         Rsq;
   Int_t           Njets_up;
   Float_t         MR_up;
   Float_t         Rsq_up;
   Int_t           Njets_down;
   Float_t         MR_down;
   Float_t         Rsq_down;
   Float_t         pileupWeight;
   Int_t           nSusyPart;
   Float_t         mSusyPart[6];   //[nSusyPart]
   Int_t           idSusyPart[6];   //[nSusyPart]
   Float_t         m22;
   Float_t         m23;
   Float_t         m24;
   Float_t         m25;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_mgg;   //!
   TBranch        *b_ptgg;   //!
   TBranch        *b_etagg;   //!
   TBranch        *b_phigg;   //!
   TBranch        *b_pho1_pt;   //!
   TBranch        *b_pho1_eta;   //!
   TBranch        *b_pho1_phi;   //!
   TBranch        *b_pho1_r9;   //!
   TBranch        *b_pho1_sigEoE;   //!
   TBranch        *b_pho1_energyGen;   //!
   TBranch        *b_pho1_genMatch;   //!
   TBranch        *b_pho1_sieie;   //!
   TBranch        *b_pho1_HE;   //!
   TBranch        *b_pho1_charged;   //!
   TBranch        *b_pho1_neutral;   //!
   TBranch        *b_pho1_photon;   //!
   TBranch        *b_pho1_eleveto;   //!
   TBranch        *b_pho1_pass_id;   //!
   TBranch        *b_pho1_pass_iso;   //!
   TBranch        *b_pho2_pt;   //!
   TBranch        *b_pho2_eta;   //!
   TBranch        *b_pho2_phi;   //!
   TBranch        *b_pho2_r9;   //!
   TBranch        *b_pho2_sigEoE;   //!
   TBranch        *b_pho2_energyGen;   //!
   TBranch        *b_pho2_genMatch;   //!
   TBranch        *b_pho2_sieie;   //!
   TBranch        *b_pho2_HE;   //!
   TBranch        *b_pho2_charged;   //!
   TBranch        *b_pho2_neutral;   //!
   TBranch        *b_pho2_photon;   //!
   TBranch        *b_pho2_eleveto;   //!
   TBranch        *b_pho2_pass_id;   //!
   TBranch        *b_pho2_pass_iso;   //!
   TBranch        *b_ele1_pt;   //!
   TBranch        *b_ele1_eta;   //!
   TBranch        *b_ele1_phi;   //!
   TBranch        *b_mu1_pt;   //!
   TBranch        *b_mu1_eta;   //!
   TBranch        *b_mu1_phi;   //!
   TBranch        *b_highest_csv;   //!
   TBranch        *b_highest_csv_pt;   //!
   TBranch        *b_highest_csv_eta;   //!
   TBranch        *b_highest_csv_phi;   //!
   TBranch        *b_highest_csv_up;   //!
   TBranch        *b_highest_csv_pt_up;   //!
   TBranch        *b_highest_csv_eta_up;   //!
   TBranch        *b_highest_csv_phi_up;   //!
   TBranch        *b_highest_csv_down;   //!
   TBranch        *b_highest_csv_pt_down;   //!
   TBranch        *b_highest_csv_eta_down;   //!
   TBranch        *b_highest_csv_phi_down;   //!
   TBranch        *b_second_csv;   //!
   TBranch        *b_second_csv_pt;   //!
   TBranch        *b_second_csv_eta;   //!
   TBranch        *b_second_csv_phi;   //!
   TBranch        *b_second_csv_up;   //!
   TBranch        *b_second_csv_pt_up;   //!
   TBranch        *b_second_csv_eta_up;   //!
   TBranch        *b_second_csv_phi_up;   //!
   TBranch        *b_second_csv_down;   //!
   TBranch        *b_second_csv_pt_down;   //!
   TBranch        *b_second_csv_eta_down;   //!
   TBranch        *b_second_csv_phi_down;   //!
   TBranch        *b_mbb;   //!
   TBranch        *b_mbb_up;   //!
   TBranch        *b_mbb_down;   //!
   TBranch        *b_mbb_NearH;   //!
   TBranch        *b_mbb_NearH_up;   //!
   TBranch        *b_mbb_NearH_down;   //!
   TBranch        *b_mbb_NearZ;   //!
   TBranch        *b_mbb_NearZ_up;   //!
   TBranch        *b_mbb_NearZ_down;   //!
   TBranch        *b_hem1_pt;   //!
   TBranch        *b_hem1_eta;   //!
   TBranch        *b_hem1_phi;   //!
   TBranch        *b_hem1_M;   //!
   TBranch        *b_hem1_pt_up;   //!
   TBranch        *b_hem1_eta_up;   //!
   TBranch        *b_hem1_phi_up;   //!
   TBranch        *b_hem1_M_up;   //!
   TBranch        *b_hem1_pt_down;   //!
   TBranch        *b_hem1_eta_down;   //!
   TBranch        *b_hem1_phi_down;   //!
   TBranch        *b_hem1_M_down;   //!
   TBranch        *b_hem2_pt;   //!
   TBranch        *b_hem2_eta;   //!
   TBranch        *b_hem2_phi;   //!
   TBranch        *b_hem2_M;   //!
   TBranch        *b_hem2_pt_up;   //!
   TBranch        *b_hem2_eta_up;   //!
   TBranch        *b_hem2_phi_up;   //!
   TBranch        *b_hem2_M_up;   //!
   TBranch        *b_hem2_pt_down;   //!
   TBranch        *b_hem2_eta_down;   //!
   TBranch        *b_hem2_phi_down;   //!
   TBranch        *b_hem2_M_down;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_MR;   //!
   TBranch        *b_Rsq;   //!
   TBranch        *b_nJets_up;   //!
   TBranch        *b_MR_up;   //!
   TBranch        *b_Rsq_up;   //!
   TBranch        *b_nJets_down;   //!
   TBranch        *b_MR_down;   //!
   TBranch        *b_Rsq_down;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_nSusyPart;   //!
   TBranch        *b_mSusyPart;   //!
   TBranch        *b_idSusyPart;   //!
   TBranch        *b_m22;   //!
   TBranch        *b_m23;   //!
   TBranch        *b_m24;   //!
   TBranch        *b_m25;   //!

   SusyHggTreeBase(TTree *tree=0);
   virtual ~SusyHggTreeBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SusyHggTreeBase_cxx
SusyHggTreeBase::SusyHggTreeBase(TTree *tree) : fChain(0) 
{
}

SusyHggTreeBase::~SusyHggTreeBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SusyHggTreeBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SusyHggTreeBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SusyHggTreeBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("mgg", &mgg, &b_mgg);
   fChain->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
   fChain->SetBranchAddress("etagg", &etagg, &b_etagg);
   fChain->SetBranchAddress("phigg", &phigg, &b_phigg);
   fChain->SetBranchAddress("pho1_pt", &pho1_pt, &b_pho1_pt);
   fChain->SetBranchAddress("pho1_eta", &pho1_eta, &b_pho1_eta);
   fChain->SetBranchAddress("pho1_phi", &pho1_phi, &b_pho1_phi);
   fChain->SetBranchAddress("pho1_r9", &pho1_r9, &b_pho1_r9);
   fChain->SetBranchAddress("pho1_sigEoE", &pho1_sigEoE, &b_pho1_sigEoE);
   fChain->SetBranchAddress("pho1_energyGen", &pho1_energyGen, &b_pho1_energyGen);
   fChain->SetBranchAddress("pho1_genMatch", &pho1_genMatch, &b_pho1_genMatch);
   fChain->SetBranchAddress("pho1_sieie", &pho1_sieie, &b_pho1_sieie);
   fChain->SetBranchAddress("pho1_HE", &pho1_HE, &b_pho1_HE);
   fChain->SetBranchAddress("pho1_charged", &pho1_charged, &b_pho1_charged);
   fChain->SetBranchAddress("pho1_neutral", &pho1_neutral, &b_pho1_neutral);
   fChain->SetBranchAddress("pho1_photon", &pho1_photon, &b_pho1_photon);
   fChain->SetBranchAddress("pho1_eleveto", &pho1_eleveto, &b_pho1_eleveto);
   fChain->SetBranchAddress("pho1_pass_id", &pho1_pass_id, &b_pho1_pass_id);
   fChain->SetBranchAddress("pho1_pass_iso", &pho1_pass_iso, &b_pho1_pass_iso);
   fChain->SetBranchAddress("pho2_pt", &pho2_pt, &b_pho2_pt);
   fChain->SetBranchAddress("pho2_eta", &pho2_eta, &b_pho2_eta);
   fChain->SetBranchAddress("pho2_phi", &pho2_phi, &b_pho2_phi);
   fChain->SetBranchAddress("pho2_r9", &pho2_r9, &b_pho2_r9);
   fChain->SetBranchAddress("pho2_sigEoE", &pho2_sigEoE, &b_pho2_sigEoE);
   fChain->SetBranchAddress("pho2_energyGen", &pho2_energyGen, &b_pho2_energyGen);
   fChain->SetBranchAddress("pho2_genMatch", &pho2_genMatch, &b_pho2_genMatch);
   fChain->SetBranchAddress("pho2_sieie", &pho2_sieie, &b_pho2_sieie);
   fChain->SetBranchAddress("pho2_HE", &pho2_HE, &b_pho2_HE);
   fChain->SetBranchAddress("pho2_charged", &pho2_charged, &b_pho2_charged);
   fChain->SetBranchAddress("pho2_neutral", &pho2_neutral, &b_pho2_neutral);
   fChain->SetBranchAddress("pho2_photon", &pho2_photon, &b_pho2_photon);
   fChain->SetBranchAddress("pho2_eleveto", &pho2_eleveto, &b_pho2_eleveto);
   fChain->SetBranchAddress("pho2_pass_id", &pho2_pass_id, &b_pho2_pass_id);
   fChain->SetBranchAddress("pho2_pass_iso", &pho2_pass_iso, &b_pho2_pass_iso);
   fChain->SetBranchAddress("ele1_pt", &ele1_pt, &b_ele1_pt);
   fChain->SetBranchAddress("ele1_eta", &ele1_eta, &b_ele1_eta);
   fChain->SetBranchAddress("ele1_phi", &ele1_phi, &b_ele1_phi);
   fChain->SetBranchAddress("mu1_pt", &mu1_pt, &b_mu1_pt);
   fChain->SetBranchAddress("mu1_eta", &mu1_eta, &b_mu1_eta);
   fChain->SetBranchAddress("mu1_phi", &mu1_phi, &b_mu1_phi);
   fChain->SetBranchAddress("highest_csv", &highest_csv, &b_highest_csv);
   fChain->SetBranchAddress("highest_csv_pt", &highest_csv_pt, &b_highest_csv_pt);
   fChain->SetBranchAddress("highest_csv_eta", &highest_csv_eta, &b_highest_csv_eta);
   fChain->SetBranchAddress("highest_csv_phi", &highest_csv_phi, &b_highest_csv_phi);
   fChain->SetBranchAddress("highest_csv_up", &highest_csv_up, &b_highest_csv_up);
   fChain->SetBranchAddress("highest_csv_pt_up", &highest_csv_pt_up, &b_highest_csv_pt_up);
   fChain->SetBranchAddress("highest_csv_eta_up", &highest_csv_eta_up, &b_highest_csv_eta_up);
   fChain->SetBranchAddress("highest_csv_phi_up", &highest_csv_phi_up, &b_highest_csv_phi_up);
   fChain->SetBranchAddress("highest_csv_down", &highest_csv_down, &b_highest_csv_down);
   fChain->SetBranchAddress("highest_csv_pt_down", &highest_csv_pt_down, &b_highest_csv_pt_down);
   fChain->SetBranchAddress("highest_csv_eta_down", &highest_csv_eta_down, &b_highest_csv_eta_down);
   fChain->SetBranchAddress("highest_csv_phi_down", &highest_csv_phi_down, &b_highest_csv_phi_down);
   fChain->SetBranchAddress("second_csv", &second_csv, &b_second_csv);
   fChain->SetBranchAddress("second_csv_pt", &second_csv_pt, &b_second_csv_pt);
   fChain->SetBranchAddress("second_csv_eta", &second_csv_eta, &b_second_csv_eta);
   fChain->SetBranchAddress("second_csv_phi", &second_csv_phi, &b_second_csv_phi);
   fChain->SetBranchAddress("second_csv_up", &second_csv_up, &b_second_csv_up);
   fChain->SetBranchAddress("second_csv_pt_up", &second_csv_pt_up, &b_second_csv_pt_up);
   fChain->SetBranchAddress("second_csv_eta_up", &second_csv_eta_up, &b_second_csv_eta_up);
   fChain->SetBranchAddress("second_csv_phi_up", &second_csv_phi_up, &b_second_csv_phi_up);
   fChain->SetBranchAddress("second_csv_down", &second_csv_down, &b_second_csv_down);
   fChain->SetBranchAddress("second_csv_pt_down", &second_csv_pt_down, &b_second_csv_pt_down);
   fChain->SetBranchAddress("second_csv_eta_down", &second_csv_eta_down, &b_second_csv_eta_down);
   fChain->SetBranchAddress("second_csv_phi_down", &second_csv_phi_down, &b_second_csv_phi_down);
   fChain->SetBranchAddress("mbb", &mbb, &b_mbb);
   fChain->SetBranchAddress("mbb_up", &mbb_up, &b_mbb_up);
   fChain->SetBranchAddress("mbb_down", &mbb_down, &b_mbb_down);
   fChain->SetBranchAddress("mbb_NearH", &mbb_NearH, &b_mbb_NearH);
   fChain->SetBranchAddress("mbb_NearH_up", &mbb_NearH_up, &b_mbb_NearH_up);
   fChain->SetBranchAddress("mbb_NearH_down", &mbb_NearH_down, &b_mbb_NearH_down);
   fChain->SetBranchAddress("mbb_NearZ", &mbb_NearZ, &b_mbb_NearZ);
   fChain->SetBranchAddress("mbb_NearZ_up", &mbb_NearZ_up, &b_mbb_NearZ_up);
   fChain->SetBranchAddress("mbb_NearZ_down", &mbb_NearZ_down, &b_mbb_NearZ_down);
   fChain->SetBranchAddress("hem1_pt", &hem1_pt, &b_hem1_pt);
   fChain->SetBranchAddress("hem1_eta", &hem1_eta, &b_hem1_eta);
   fChain->SetBranchAddress("hem1_phi", &hem1_phi, &b_hem1_phi);
   fChain->SetBranchAddress("hem1_M", &hem1_M, &b_hem1_M);
   fChain->SetBranchAddress("hem1_pt_up", &hem1_pt_up, &b_hem1_pt_up);
   fChain->SetBranchAddress("hem1_eta_up", &hem1_eta_up, &b_hem1_eta_up);
   fChain->SetBranchAddress("hem1_phi_up", &hem1_phi_up, &b_hem1_phi_up);
   fChain->SetBranchAddress("hem1_M_up", &hem1_M_up, &b_hem1_M_up);
   fChain->SetBranchAddress("hem1_pt_down", &hem1_pt_down, &b_hem1_pt_down);
   fChain->SetBranchAddress("hem1_eta_down", &hem1_eta_down, &b_hem1_eta_down);
   fChain->SetBranchAddress("hem1_phi_down", &hem1_phi_down, &b_hem1_phi_down);
   fChain->SetBranchAddress("hem1_M_down", &hem1_M_down, &b_hem1_M_down);
   fChain->SetBranchAddress("hem2_pt", &hem2_pt, &b_hem2_pt);
   fChain->SetBranchAddress("hem2_eta", &hem2_eta, &b_hem2_eta);
   fChain->SetBranchAddress("hem2_phi", &hem2_phi, &b_hem2_phi);
   fChain->SetBranchAddress("hem2_M", &hem2_M, &b_hem2_M);
   fChain->SetBranchAddress("hem2_pt_up", &hem2_pt_up, &b_hem2_pt_up);
   fChain->SetBranchAddress("hem2_eta_up", &hem2_eta_up, &b_hem2_eta_up);
   fChain->SetBranchAddress("hem2_phi_up", &hem2_phi_up, &b_hem2_phi_up);
   fChain->SetBranchAddress("hem2_M_up", &hem2_M_up, &b_hem2_M_up);
   fChain->SetBranchAddress("hem2_pt_down", &hem2_pt_down, &b_hem2_pt_down);
   fChain->SetBranchAddress("hem2_eta_down", &hem2_eta_down, &b_hem2_eta_down);
   fChain->SetBranchAddress("hem2_phi_down", &hem2_phi_down, &b_hem2_phi_down);
   fChain->SetBranchAddress("hem2_M_down", &hem2_M_down, &b_hem2_M_down);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("Njets", &Njets, &b_nJets);
   fChain->SetBranchAddress("MR", &MR, &b_MR);
   fChain->SetBranchAddress("Rsq", &Rsq, &b_Rsq);
   fChain->SetBranchAddress("Njets_up", &Njets_up, &b_nJets_up);
   fChain->SetBranchAddress("MR_up", &MR_up, &b_MR_up);
   fChain->SetBranchAddress("Rsq_up", &Rsq_up, &b_Rsq_up);
   fChain->SetBranchAddress("Njets_down", &Njets_down, &b_nJets_down);
   fChain->SetBranchAddress("MR_down", &MR_down, &b_MR_down);
   fChain->SetBranchAddress("Rsq_down", &Rsq_down, &b_Rsq_down);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("nSusyPart", &nSusyPart, &b_nSusyPart);
   fChain->SetBranchAddress("mSusyPart", mSusyPart, &b_mSusyPart);
   fChain->SetBranchAddress("idSusyPart", idSusyPart, &b_idSusyPart);
   fChain->SetBranchAddress("m22", &m22, &b_m22);
   fChain->SetBranchAddress("m23", &m23, &b_m23);
   fChain->SetBranchAddress("m24", &m24, &b_m24);
   fChain->SetBranchAddress("m25", &m25, &b_m25);
   Notify();
}

Bool_t SusyHggTreeBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SusyHggTreeBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SusyHggTreeBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SusyHggTreeBase_cxx
