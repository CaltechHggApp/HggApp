//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov  5 00:33:04 2012 by ROOT version 5.32/00
// from TTree HggOutput/
// found on file: DiPhotonJets_madgraph_DR13X-PU_S10_START53.root
//////////////////////////////////////////////////////////

#ifndef HggOutputReader_h
#define HggOutputReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.
#ifndef HggOutputReader2_h
const Int_t kMaxPhoton = 2;
const Int_t kMaxPhotonPFCiC = 2;
const Int_t kMaxPhotonCiC = 2;
#endif

class HggOutputReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           trigger;
   Float_t         mPair;
   Float_t         mPairNoCorr;
   Float_t         mPairRes;
   Float_t         mPairResWrongVtx;
   Float_t         diPhotonMVA;
   Int_t           diPhotonVtx;
   Float_t         diPhotonVtxX;
   Float_t         diPhotonVtxY;
   Float_t         diPhotonVtxZ;
   Float_t         vtxProb;
   Float_t         Mjj;
   Float_t         ptJet1;
   Float_t         ptJet2;
   Float_t         cosThetaLead;
   Int_t           cat;
   Float_t         mPairPFCiC;
   Float_t         mPairNoCorrPFCiC;
   Float_t         mPairResPFCiC;
   Float_t         mPairResWrongVtxPFCiC;
   Int_t           diPhotonVtxPFCiC;
   Float_t         diPhotonVtxXPFCiC;
   Float_t         diPhotonVtxYPFCiC;
   Float_t         diPhotonVtxZPFCiC;
   Float_t         vtxProbPFCiC;
   Float_t         MjjPFCiC;
   Float_t         ptJet1PFCiC;
   Float_t         ptJet2PFCiC;
   Float_t         cosThetaLeadPFCiC;
   Int_t           catPFCiC;
   Float_t         mPairCiC;
   Float_t         mPairNoCorrCiC;
   Float_t         mPairResCiC;
   Float_t         mPairResWrongVtxCiC;
   Int_t           diPhotonVtxCiC;
   Float_t         diPhotonVtxXCiC;
   Float_t         diPhotonVtxYCiC;
   Float_t         diPhotonVtxZCiC;
   Float_t         vtxProbCiC;
   Float_t         MjjCiC;
   Float_t         ptJet1CiC;
   Float_t         ptJet2CiC;
   Float_t         cosThetaLeadCiC;
   Int_t           nPhoton;
   Int_t           Photon_;
   Float_t         Photon_pt[kMaxPhoton];   //[Photon_]
   Float_t         Photon_eta[kMaxPhoton];   //[Photon_]
   Float_t         Photon_phi[kMaxPhoton];   //[Photon_]
   Float_t         Photon_E[kMaxPhoton];   //[Photon_]
   Float_t         Photon_EError[kMaxPhoton];   //[Photon_]
   Float_t         Photon_EErrorSmeared[kMaxPhoton];   //[Photon_]
   Float_t         Photon_pt_Gen[kMaxPhoton];   //[Photon_]
   Float_t         Photon_eta_Gen[kMaxPhoton];   //[Photon_]
   Float_t         Photon_phi_Gen[kMaxPhoton];   //[Photon_]
   Float_t         Photon_E_Gen[kMaxPhoton];   //[Photon_]
   Float_t         Photon_etaSC[kMaxPhoton];   //[Photon_]
   Int_t           Photon_index[kMaxPhoton];   //[Photon_]
   Float_t         Photon_r9[kMaxPhoton];   //[Photon_]
   Bool_t          Photon_passPFCiC[kMaxPhoton];   //[Photon_]
   Int_t           Photon_category[kMaxPhoton];   //[Photon_]
   Float_t         Photon_idMVA[kMaxPhoton];   //[Photon_]
   Float_t         Photon_HoverE[kMaxPhoton];   //[Photon_]
   Float_t         Photon_sieie[kMaxPhoton];   //[Photon_]
   Float_t         Photon_dr03PFChargedIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_isosumGood[kMaxPhoton];   //[Photon_]
   Float_t         Photon_isosumBad[kMaxPhoton];   //[Photon_]
   Float_t         Photon_dr03EcalIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_dr04HcalIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_dr03TrackIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_dr02PFChargedIso[kMaxPhoton];   //[Photon_]
   Int_t           nPhotonPFCiC;
   Int_t           PhotonPFCiC_;
   Float_t         PhotonPFCiC_pt[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_eta[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_phi[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_E[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_EError[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_EErrorSmeared[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_pt_Gen[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_eta_Gen[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_phi_Gen[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_E_Gen[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_etaSC[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Int_t           PhotonPFCiC_index[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_r9[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Bool_t          PhotonPFCiC_passPFCiC[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Int_t           PhotonPFCiC_category[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_idMVA[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_HoverE[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_sieie[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_dr03PFChargedIso[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_isosumGood[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_isosumBad[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_dr03EcalIso[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_dr04HcalIso[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_dr03TrackIso[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Float_t         PhotonPFCiC_dr02PFChargedIso[kMaxPhotonPFCiC];   //[PhotonPFCiC_]
   Int_t           nPhotonCiC;
   Int_t           PhotonCiC_;
   Float_t         PhotonCiC_pt[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_eta[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_phi[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_E[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_EError[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_EErrorSmeared[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_pt_Gen[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_eta_Gen[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_phi_Gen[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_E_Gen[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_etaSC[kMaxPhotonCiC];   //[PhotonCiC_]
   Int_t           PhotonCiC_index[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_r9[kMaxPhotonCiC];   //[PhotonCiC_]
   Bool_t          PhotonCiC_passPFCiC[kMaxPhotonCiC];   //[PhotonCiC_]
   Int_t           PhotonCiC_category[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_idMVA[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_HoverE[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_sieie[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_dr03PFChargedIso[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_isosumGood[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_isosumBad[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_dr03EcalIso[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_dr04HcalIso[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_dr03TrackIso[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         PhotonCiC_dr02PFChargedIso[kMaxPhotonCiC];   //[PhotonCiC_]
   Float_t         MET;
   Float_t         METPhi;
   Float_t         genHiggsPt;
   Float_t         genHiggsVx;
   Float_t         genHiggsVy;
   Float_t         genHiggsVz;
   Float_t         evtWeight;
   Float_t         ptGenPho1;
   Float_t         etaGenPho1;
   Float_t         phiGenPho1;
   Float_t         energyGenPho1;
   Float_t         ptGenPho2;
   Float_t         etaGenPho2;
   Float_t         phiGenPho2;
   Float_t         energyGenPho2;
   vector<float>   *mPairScale;
   vector<float>   *pho1MVAScale;
   vector<float>   *pho2MVAScale;
   vector<float>   *diPhoMVAScale;
   vector<float>   *mPairSmear;
   vector<float>   *pho1MVASmear;
   vector<float>   *pho2MVASmear;
   vector<float>   *diPhoMVASmear;
   vector<float>   *mPairScalePFCiC;
   vector<float>   *mPairSmearPFCiC;
   vector<float>   *mPairScaleCiC;
   vector<float>   *mPairSmearCiC;
   Float_t         nPU;
   Int_t           nVtx;
   Int_t           runNumber;
   Int_t           evtNumber;
   Int_t           lumiBlock;

   // List of branches
   TBranch        *b_trigger;   //!
   TBranch        *b_mPair;   //!
   TBranch        *b_mPairNoCorr;   //!
   TBranch        *b_mPairRes;   //!
   TBranch        *b_mPairResWrongVtx;   //!
   TBranch        *b_diPhotonMVA;   //!
   TBranch        *b_diPhotonVtx;   //!
   TBranch        *b_diPhotonVtxX;   //!
   TBranch        *b_diPhotonVtxY;   //!
   TBranch        *b_diPhotonVtxZ;   //!
   TBranch        *b_vtxProb;   //!
   TBranch        *b_Mjj;   //!
   TBranch        *b_ptJet1;   //!
   TBranch        *b_ptJet2;   //!
   TBranch        *b_cosThetaLead;   //!
   TBranch        *b_cat;   //!
   TBranch        *b_mPairPFCiC;   //!
   TBranch        *b_mPairNoCorrPFCiC;   //!
   TBranch        *b_mPairResPFCiC;   //!
   TBranch        *b_mPairResWrongVtxPFCiC;   //!
   TBranch        *b_diPhotonVtxPFCiC;   //!
   TBranch        *b_diPhotonVtxXPFCiC;   //!
   TBranch        *b_diPhotonVtxYPFCiC;   //!
   TBranch        *b_diPhotonVtxZPFCiC;   //!
   TBranch        *b_vtxProbPFCiC;   //!
   TBranch        *b_MjjPFCiC;   //!
   TBranch        *b_ptJet1PFCiC;   //!
   TBranch        *b_ptJet2PFCiC;   //!
   TBranch        *b_cosThetaLeadPFCiC;   //!
   TBranch        *b_catPFCiC;   //!
   TBranch        *b_mPairCiC;   //!
   TBranch        *b_mPairNoCorrCiC;   //!
   TBranch        *b_mPairResCiC;   //!
   TBranch        *b_mPairResWrongVtxCiC;   //!
   TBranch        *b_diPhotonVtxCiC;   //!
   TBranch        *b_diPhotonVtxXCiC;   //!
   TBranch        *b_diPhotonVtxYCiC;   //!
   TBranch        *b_diPhotonVtxZCiC;   //!
   TBranch        *b_vtxProbCiC;   //!
   TBranch        *b_MjjCiC;   //!
   TBranch        *b_ptJet1CiC;   //!
   TBranch        *b_ptJet2CiC;   //!
   TBranch        *b_cosThetaLeadCiC;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_EError;   //!
   TBranch        *b_Photon_EErrorSmeared;   //!
   TBranch        *b_Photon_pt_Gen;   //!
   TBranch        *b_Photon_eta_Gen;   //!
   TBranch        *b_Photon_phi_Gen;   //!
   TBranch        *b_Photon_E_Gen;   //!
   TBranch        *b_Photon_etaSC;   //!
   TBranch        *b_Photon_index;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_passPFCiC;   //!
   TBranch        *b_Photon_category;   //!
   TBranch        *b_Photon_idMVA;   //!
   TBranch        *b_Photon_HoverE;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_dr03PFChargedIso;   //!
   TBranch        *b_Photon_isosumGood;   //!
   TBranch        *b_Photon_isosumBad;   //!
   TBranch        *b_Photon_dr03EcalIso;   //!
   TBranch        *b_Photon_dr04HcalIso;   //!
   TBranch        *b_Photon_dr03TrackIso;   //!
   TBranch        *b_Photon_dr02PFChargedIso;   //!
   TBranch        *b_nPhotonPFCiC;   //!
   TBranch        *b_PhotonPFCiC_;   //!
   TBranch        *b_PhotonPFCiC_pt;   //!
   TBranch        *b_PhotonPFCiC_eta;   //!
   TBranch        *b_PhotonPFCiC_phi;   //!
   TBranch        *b_PhotonPFCiC_E;   //!
   TBranch        *b_PhotonPFCiC_EError;   //!
   TBranch        *b_PhotonPFCiC_EErrorSmeared;   //!
   TBranch        *b_PhotonPFCiC_pt_Gen;   //!
   TBranch        *b_PhotonPFCiC_eta_Gen;   //!
   TBranch        *b_PhotonPFCiC_phi_Gen;   //!
   TBranch        *b_PhotonPFCiC_E_Gen;   //!
   TBranch        *b_PhotonPFCiC_etaSC;   //!
   TBranch        *b_PhotonPFCiC_index;   //!
   TBranch        *b_PhotonPFCiC_r9;   //!
   TBranch        *b_PhotonPFCiC_passPFCiC;   //!
   TBranch        *b_PhotonPFCiC_category;   //!
   TBranch        *b_PhotonPFCiC_idMVA;   //!
   TBranch        *b_PhotonPFCiC_HoverE;   //!
   TBranch        *b_PhotonPFCiC_sieie;   //!
   TBranch        *b_PhotonPFCiC_dr03PFChargedIso;   //!
   TBranch        *b_PhotonPFCiC_isosumGood;   //!
   TBranch        *b_PhotonPFCiC_isosumBad;   //!
   TBranch        *b_PhotonPFCiC_dr03EcalIso;   //!
   TBranch        *b_PhotonPFCiC_dr04HcalIso;   //!
   TBranch        *b_PhotonPFCiC_dr03TrackIso;   //!
   TBranch        *b_PhotonPFCiC_dr02PFChargedIso;   //!
   TBranch        *b_nPhotonCiC;   //!
   TBranch        *b_PhotonCiC_;   //!
   TBranch        *b_PhotonCiC_pt;   //!
   TBranch        *b_PhotonCiC_eta;   //!
   TBranch        *b_PhotonCiC_phi;   //!
   TBranch        *b_PhotonCiC_E;   //!
   TBranch        *b_PhotonCiC_EError;   //!
   TBranch        *b_PhotonCiC_EErrorSmeared;   //!
   TBranch        *b_PhotonCiC_pt_Gen;   //!
   TBranch        *b_PhotonCiC_eta_Gen;   //!
   TBranch        *b_PhotonCiC_phi_Gen;   //!
   TBranch        *b_PhotonCiC_E_Gen;   //!
   TBranch        *b_PhotonCiC_etaSC;   //!
   TBranch        *b_PhotonCiC_index;   //!
   TBranch        *b_PhotonCiC_r9;   //!
   TBranch        *b_PhotonCiC_passPFCiC;   //!
   TBranch        *b_PhotonCiC_category;   //!
   TBranch        *b_PhotonCiC_idMVA;   //!
   TBranch        *b_PhotonCiC_HoverE;   //!
   TBranch        *b_PhotonCiC_sieie;   //!
   TBranch        *b_PhotonCiC_dr03PFChargedIso;   //!
   TBranch        *b_PhotonCiC_isosumGood;   //!
   TBranch        *b_PhotonCiC_isosumBad;   //!
   TBranch        *b_PhotonCiC_dr03EcalIso;   //!
   TBranch        *b_PhotonCiC_dr04HcalIso;   //!
   TBranch        *b_PhotonCiC_dr03TrackIso;   //!
   TBranch        *b_PhotonCiC_dr02PFChargedIso;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_genHiggsPt;   //!
   TBranch        *b_genHiggsVx;   //!
   TBranch        *b_genHiggsVy;   //!
   TBranch        *b_genHiggsVz;   //!
   TBranch        *b_evtWeight;   //!
   TBranch        *b_ptGenPho1;   //!
   TBranch        *b_etaGenPho1;   //!
   TBranch        *b_phiGenPho1;   //!
   TBranch        *b_energyGenPho1;   //!
   TBranch        *b_ptGenPho2;   //!
   TBranch        *b_etaGenPho2;   //!
   TBranch        *b_phiGenPho2;   //!
   TBranch        *b_energyGenPho2;   //!
   TBranch        *b_mPairScale;   //!
   TBranch        *b_pho1MVAScale;   //!
   TBranch        *b_pho2MVAScale;   //!
   TBranch        *b_diPhoMVAScale;   //!
   TBranch        *b_mPairSmear;   //!
   TBranch        *b_pho1MVASmear;   //!
   TBranch        *b_pho2MVASmear;   //!
   TBranch        *b_diPhoMVASmear;   //!
   TBranch        *b_mPairScalePFCiC;   //!
   TBranch        *b_mPairSmearPFCiC;   //!
   TBranch        *b_mPairScaleCiC;   //!
   TBranch        *b_mPairSmearCiC;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_evtNumber;   //!
   TBranch        *b_lumiBlock;   //!

   HggOutputReader(TTree *tree=0);
   virtual ~HggOutputReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HggOutputReader_cxx
HggOutputReader::HggOutputReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DiPhotonJets_madgraph_DR13X-PU_S10_START53.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DiPhotonJets_madgraph_DR13X-PU_S10_START53.root");
      }
      f->GetObject("HggOutput",tree);

   }
   Init(tree);
}

HggOutputReader::~HggOutputReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HggOutputReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HggOutputReader::LoadTree(Long64_t entry)
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

void HggOutputReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mPairScale = 0;
   pho1MVAScale = 0;
   pho2MVAScale = 0;
   diPhoMVAScale = 0;
   mPairSmear = 0;
   pho1MVASmear = 0;
   pho2MVASmear = 0;
   diPhoMVASmear = 0;
   mPairScalePFCiC = 0;
   mPairSmearPFCiC = 0;
   mPairScaleCiC = 0;
   mPairSmearCiC = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("mPair", &mPair, &b_mPair);
   fChain->SetBranchAddress("mPairNoCorr", &mPairNoCorr, &b_mPairNoCorr);
   fChain->SetBranchAddress("mPairRes", &mPairRes, &b_mPairRes);
   fChain->SetBranchAddress("mPairResWrongVtx", &mPairResWrongVtx, &b_mPairResWrongVtx);
   fChain->SetBranchAddress("diPhotonMVA", &diPhotonMVA, &b_diPhotonMVA);
   fChain->SetBranchAddress("diPhotonVtx", &diPhotonVtx, &b_diPhotonVtx);
   fChain->SetBranchAddress("diPhotonVtxX", &diPhotonVtxX, &b_diPhotonVtxX);
   fChain->SetBranchAddress("diPhotonVtxY", &diPhotonVtxY, &b_diPhotonVtxY);
   fChain->SetBranchAddress("diPhotonVtxZ", &diPhotonVtxZ, &b_diPhotonVtxZ);
   fChain->SetBranchAddress("vtxProb", &vtxProb, &b_vtxProb);
   fChain->SetBranchAddress("Mjj", &Mjj, &b_Mjj);
   fChain->SetBranchAddress("ptJet1", &ptJet1, &b_ptJet1);
   fChain->SetBranchAddress("ptJet2", &ptJet2, &b_ptJet2);
   fChain->SetBranchAddress("cosThetaLead", &cosThetaLead, &b_cosThetaLead);
   fChain->SetBranchAddress("cat", &cat, &b_cat);
   fChain->SetBranchAddress("mPairPFCiC", &mPairPFCiC, &b_mPairPFCiC);
   fChain->SetBranchAddress("mPairNoCorrPFCiC", &mPairNoCorrPFCiC, &b_mPairNoCorrPFCiC);
   fChain->SetBranchAddress("mPairResPFCiC", &mPairResPFCiC, &b_mPairResPFCiC);
   fChain->SetBranchAddress("mPairResWrongVtxPFCiC", &mPairResWrongVtxPFCiC, &b_mPairResWrongVtxPFCiC);
   fChain->SetBranchAddress("diPhotonVtxPFCiC", &diPhotonVtxPFCiC, &b_diPhotonVtxPFCiC);
   fChain->SetBranchAddress("diPhotonVtxXPFCiC", &diPhotonVtxXPFCiC, &b_diPhotonVtxXPFCiC);
   fChain->SetBranchAddress("diPhotonVtxYPFCiC", &diPhotonVtxYPFCiC, &b_diPhotonVtxYPFCiC);
   fChain->SetBranchAddress("diPhotonVtxZPFCiC", &diPhotonVtxZPFCiC, &b_diPhotonVtxZPFCiC);
   fChain->SetBranchAddress("vtxProbPFCiC", &vtxProbPFCiC, &b_vtxProbPFCiC);
   fChain->SetBranchAddress("MjjPFCiC", &MjjPFCiC, &b_MjjPFCiC);
   fChain->SetBranchAddress("ptJet1PFCiC", &ptJet1PFCiC, &b_ptJet1PFCiC);
   fChain->SetBranchAddress("ptJet2PFCiC", &ptJet2PFCiC, &b_ptJet2PFCiC);
   fChain->SetBranchAddress("cosThetaLeadPFCiC", &cosThetaLeadPFCiC, &b_cosThetaLeadPFCiC);
   fChain->SetBranchAddress("catPFCiC", &catPFCiC, &b_catPFCiC);
   fChain->SetBranchAddress("mPairCiC", &mPairCiC, &b_mPairCiC);
   fChain->SetBranchAddress("mPairNoCorrCiC", &mPairNoCorrCiC, &b_mPairNoCorrCiC);
   fChain->SetBranchAddress("mPairResCiC", &mPairResCiC, &b_mPairResCiC);
   fChain->SetBranchAddress("mPairResWrongVtxCiC", &mPairResWrongVtxCiC, &b_mPairResWrongVtxCiC);
   fChain->SetBranchAddress("diPhotonVtxCiC", &diPhotonVtxCiC, &b_diPhotonVtxCiC);
   fChain->SetBranchAddress("diPhotonVtxXCiC", &diPhotonVtxXCiC, &b_diPhotonVtxXCiC);
   fChain->SetBranchAddress("diPhotonVtxYCiC", &diPhotonVtxYCiC, &b_diPhotonVtxYCiC);
   fChain->SetBranchAddress("diPhotonVtxZCiC", &diPhotonVtxZCiC, &b_diPhotonVtxZCiC);
   fChain->SetBranchAddress("vtxProbCiC", &vtxProbCiC, &b_vtxProbCiC);
   fChain->SetBranchAddress("MjjCiC", &MjjCiC, &b_MjjCiC);
   fChain->SetBranchAddress("ptJet1CiC", &ptJet1CiC, &b_ptJet1CiC);
   fChain->SetBranchAddress("ptJet2CiC", &ptJet2CiC, &b_ptJet2CiC);
   fChain->SetBranchAddress("cosThetaLeadCiC", &cosThetaLeadCiC, &b_cosThetaLeadCiC);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon", &Photon_, &b_Photon_);
   fChain->SetBranchAddress("Photon.pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon.eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon.phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon.E", Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon.EError", Photon_EError, &b_Photon_EError);
   fChain->SetBranchAddress("Photon.EErrorSmeared", Photon_EErrorSmeared, &b_Photon_EErrorSmeared);
   fChain->SetBranchAddress("Photon.pt_Gen", Photon_pt_Gen, &b_Photon_pt_Gen);
   fChain->SetBranchAddress("Photon.eta_Gen", Photon_eta_Gen, &b_Photon_eta_Gen);
   fChain->SetBranchAddress("Photon.phi_Gen", Photon_phi_Gen, &b_Photon_phi_Gen);
   fChain->SetBranchAddress("Photon.E_Gen", Photon_E_Gen, &b_Photon_E_Gen);
   fChain->SetBranchAddress("Photon.etaSC", Photon_etaSC, &b_Photon_etaSC);
   fChain->SetBranchAddress("Photon.index", Photon_index, &b_Photon_index);
   fChain->SetBranchAddress("Photon.r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon.passPFCiC", Photon_passPFCiC, &b_Photon_passPFCiC);
   fChain->SetBranchAddress("Photon.category", Photon_category, &b_Photon_category);
   fChain->SetBranchAddress("Photon.idMVA", Photon_idMVA, &b_Photon_idMVA);
   fChain->SetBranchAddress("Photon.HoverE", Photon_HoverE, &b_Photon_HoverE);
   fChain->SetBranchAddress("Photon.sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon.dr03PFChargedIso", Photon_dr03PFChargedIso, &b_Photon_dr03PFChargedIso);
   fChain->SetBranchAddress("Photon.isosumGood", Photon_isosumGood, &b_Photon_isosumGood);
   fChain->SetBranchAddress("Photon.isosumBad", Photon_isosumBad, &b_Photon_isosumBad);
   fChain->SetBranchAddress("Photon.dr03EcalIso", Photon_dr03EcalIso, &b_Photon_dr03EcalIso);
   fChain->SetBranchAddress("Photon.dr04HcalIso", Photon_dr04HcalIso, &b_Photon_dr04HcalIso);
   fChain->SetBranchAddress("Photon.dr03TrackIso", Photon_dr03TrackIso, &b_Photon_dr03TrackIso);
   fChain->SetBranchAddress("Photon.dr02PFChargedIso", Photon_dr02PFChargedIso, &b_Photon_dr02PFChargedIso);
   fChain->SetBranchAddress("nPhotonPFCiC", &nPhotonPFCiC, &b_nPhotonPFCiC);
   fChain->SetBranchAddress("PhotonPFCiC", &PhotonPFCiC_, &b_PhotonPFCiC_);
   fChain->SetBranchAddress("PhotonPFCiC.pt", PhotonPFCiC_pt, &b_PhotonPFCiC_pt);
   fChain->SetBranchAddress("PhotonPFCiC.eta", PhotonPFCiC_eta, &b_PhotonPFCiC_eta);
   fChain->SetBranchAddress("PhotonPFCiC.phi", PhotonPFCiC_phi, &b_PhotonPFCiC_phi);
   fChain->SetBranchAddress("PhotonPFCiC.E", PhotonPFCiC_E, &b_PhotonPFCiC_E);
   fChain->SetBranchAddress("PhotonPFCiC.EError", PhotonPFCiC_EError, &b_PhotonPFCiC_EError);
   fChain->SetBranchAddress("PhotonPFCiC.EErrorSmeared", PhotonPFCiC_EErrorSmeared, &b_PhotonPFCiC_EErrorSmeared);
   fChain->SetBranchAddress("PhotonPFCiC.pt_Gen", PhotonPFCiC_pt_Gen, &b_PhotonPFCiC_pt_Gen);
   fChain->SetBranchAddress("PhotonPFCiC.eta_Gen", PhotonPFCiC_eta_Gen, &b_PhotonPFCiC_eta_Gen);
   fChain->SetBranchAddress("PhotonPFCiC.phi_Gen", PhotonPFCiC_phi_Gen, &b_PhotonPFCiC_phi_Gen);
   fChain->SetBranchAddress("PhotonPFCiC.E_Gen", PhotonPFCiC_E_Gen, &b_PhotonPFCiC_E_Gen);
   fChain->SetBranchAddress("PhotonPFCiC.etaSC", PhotonPFCiC_etaSC, &b_PhotonPFCiC_etaSC);
   fChain->SetBranchAddress("PhotonPFCiC.index", PhotonPFCiC_index, &b_PhotonPFCiC_index);
   fChain->SetBranchAddress("PhotonPFCiC.r9", PhotonPFCiC_r9, &b_PhotonPFCiC_r9);
   fChain->SetBranchAddress("PhotonPFCiC.passPFCiC", PhotonPFCiC_passPFCiC, &b_PhotonPFCiC_passPFCiC);
   fChain->SetBranchAddress("PhotonPFCiC.category", PhotonPFCiC_category, &b_PhotonPFCiC_category);
   fChain->SetBranchAddress("PhotonPFCiC.idMVA", PhotonPFCiC_idMVA, &b_PhotonPFCiC_idMVA);
   fChain->SetBranchAddress("PhotonPFCiC.HoverE", PhotonPFCiC_HoverE, &b_PhotonPFCiC_HoverE);
   fChain->SetBranchAddress("PhotonPFCiC.sieie", PhotonPFCiC_sieie, &b_PhotonPFCiC_sieie);
   fChain->SetBranchAddress("PhotonPFCiC.dr03PFChargedIso", PhotonPFCiC_dr03PFChargedIso, &b_PhotonPFCiC_dr03PFChargedIso);
   fChain->SetBranchAddress("PhotonPFCiC.isosumGood", PhotonPFCiC_isosumGood, &b_PhotonPFCiC_isosumGood);
   fChain->SetBranchAddress("PhotonPFCiC.isosumBad", PhotonPFCiC_isosumBad, &b_PhotonPFCiC_isosumBad);
   fChain->SetBranchAddress("PhotonPFCiC.dr03EcalIso", PhotonPFCiC_dr03EcalIso, &b_PhotonPFCiC_dr03EcalIso);
   fChain->SetBranchAddress("PhotonPFCiC.dr04HcalIso", PhotonPFCiC_dr04HcalIso, &b_PhotonPFCiC_dr04HcalIso);
   fChain->SetBranchAddress("PhotonPFCiC.dr03TrackIso", PhotonPFCiC_dr03TrackIso, &b_PhotonPFCiC_dr03TrackIso);
   fChain->SetBranchAddress("PhotonPFCiC.dr02PFChargedIso", PhotonPFCiC_dr02PFChargedIso, &b_PhotonPFCiC_dr02PFChargedIso);
   fChain->SetBranchAddress("nPhotonCiC", &nPhotonCiC, &b_nPhotonCiC);
   fChain->SetBranchAddress("PhotonCiC", &PhotonCiC_, &b_PhotonCiC_);
   fChain->SetBranchAddress("PhotonCiC.pt", PhotonCiC_pt, &b_PhotonCiC_pt);
   fChain->SetBranchAddress("PhotonCiC.eta", PhotonCiC_eta, &b_PhotonCiC_eta);
   fChain->SetBranchAddress("PhotonCiC.phi", PhotonCiC_phi, &b_PhotonCiC_phi);
   fChain->SetBranchAddress("PhotonCiC.E", PhotonCiC_E, &b_PhotonCiC_E);
   fChain->SetBranchAddress("PhotonCiC.EError", PhotonCiC_EError, &b_PhotonCiC_EError);
   fChain->SetBranchAddress("PhotonCiC.EErrorSmeared", PhotonCiC_EErrorSmeared, &b_PhotonCiC_EErrorSmeared);
   fChain->SetBranchAddress("PhotonCiC.pt_Gen", PhotonCiC_pt_Gen, &b_PhotonCiC_pt_Gen);
   fChain->SetBranchAddress("PhotonCiC.eta_Gen", PhotonCiC_eta_Gen, &b_PhotonCiC_eta_Gen);
   fChain->SetBranchAddress("PhotonCiC.phi_Gen", PhotonCiC_phi_Gen, &b_PhotonCiC_phi_Gen);
   fChain->SetBranchAddress("PhotonCiC.E_Gen", PhotonCiC_E_Gen, &b_PhotonCiC_E_Gen);
   fChain->SetBranchAddress("PhotonCiC.etaSC", PhotonCiC_etaSC, &b_PhotonCiC_etaSC);
   fChain->SetBranchAddress("PhotonCiC.index", PhotonCiC_index, &b_PhotonCiC_index);
   fChain->SetBranchAddress("PhotonCiC.r9", PhotonCiC_r9, &b_PhotonCiC_r9);
   fChain->SetBranchAddress("PhotonCiC.passPFCiC", PhotonCiC_passPFCiC, &b_PhotonCiC_passPFCiC);
   fChain->SetBranchAddress("PhotonCiC.category", PhotonCiC_category, &b_PhotonCiC_category);
   fChain->SetBranchAddress("PhotonCiC.idMVA", PhotonCiC_idMVA, &b_PhotonCiC_idMVA);
   fChain->SetBranchAddress("PhotonCiC.HoverE", PhotonCiC_HoverE, &b_PhotonCiC_HoverE);
   fChain->SetBranchAddress("PhotonCiC.sieie", PhotonCiC_sieie, &b_PhotonCiC_sieie);
   fChain->SetBranchAddress("PhotonCiC.dr03PFChargedIso", PhotonCiC_dr03PFChargedIso, &b_PhotonCiC_dr03PFChargedIso);
   fChain->SetBranchAddress("PhotonCiC.isosumGood", PhotonCiC_isosumGood, &b_PhotonCiC_isosumGood);
   fChain->SetBranchAddress("PhotonCiC.isosumBad", PhotonCiC_isosumBad, &b_PhotonCiC_isosumBad);
   fChain->SetBranchAddress("PhotonCiC.dr03EcalIso", PhotonCiC_dr03EcalIso, &b_PhotonCiC_dr03EcalIso);
   fChain->SetBranchAddress("PhotonCiC.dr04HcalIso", PhotonCiC_dr04HcalIso, &b_PhotonCiC_dr04HcalIso);
   fChain->SetBranchAddress("PhotonCiC.dr03TrackIso", PhotonCiC_dr03TrackIso, &b_PhotonCiC_dr03TrackIso);
   fChain->SetBranchAddress("PhotonCiC.dr02PFChargedIso", PhotonCiC_dr02PFChargedIso, &b_PhotonCiC_dr02PFChargedIso);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("genHiggsPt", &genHiggsPt, &b_genHiggsPt);
   fChain->SetBranchAddress("genHiggsVx", &genHiggsVx, &b_genHiggsVx);
   fChain->SetBranchAddress("genHiggsVy", &genHiggsVy, &b_genHiggsVy);
   fChain->SetBranchAddress("genHiggsVz", &genHiggsVz, &b_genHiggsVz);
   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
   fChain->SetBranchAddress("ptGenPho1", &ptGenPho1, &b_ptGenPho1);
   fChain->SetBranchAddress("etaGenPho1", &etaGenPho1, &b_etaGenPho1);
   fChain->SetBranchAddress("phiGenPho1", &phiGenPho1, &b_phiGenPho1);
   fChain->SetBranchAddress("energyGenPho1", &energyGenPho1, &b_energyGenPho1);
   fChain->SetBranchAddress("ptGenPho2", &ptGenPho2, &b_ptGenPho2);
   fChain->SetBranchAddress("etaGenPho2", &etaGenPho2, &b_etaGenPho2);
   fChain->SetBranchAddress("phiGenPho2", &phiGenPho2, &b_phiGenPho2);
   fChain->SetBranchAddress("energyGenPho2", &energyGenPho2, &b_energyGenPho2);
   fChain->SetBranchAddress("mPairScale", &mPairScale, &b_mPairScale);
   fChain->SetBranchAddress("pho1MVAScale", &pho1MVAScale, &b_pho1MVAScale);
   fChain->SetBranchAddress("pho2MVAScale", &pho2MVAScale, &b_pho2MVAScale);
   fChain->SetBranchAddress("diPhoMVAScale", &diPhoMVAScale, &b_diPhoMVAScale);
   fChain->SetBranchAddress("mPairSmear", &mPairSmear, &b_mPairSmear);
   fChain->SetBranchAddress("pho1MVASmear", &pho1MVASmear, &b_pho1MVASmear);
   fChain->SetBranchAddress("pho2MVASmear", &pho2MVASmear, &b_pho2MVASmear);
   fChain->SetBranchAddress("diPhoMVASmear", &diPhoMVASmear, &b_diPhoMVASmear);
   fChain->SetBranchAddress("mPairScalePFCiC", &mPairScalePFCiC, &b_mPairScalePFCiC);
   fChain->SetBranchAddress("mPairSmearPFCiC", &mPairSmearPFCiC, &b_mPairSmearPFCiC);
   fChain->SetBranchAddress("mPairScaleCiC", &mPairScaleCiC, &b_mPairScaleCiC);
   fChain->SetBranchAddress("mPairSmearCiC", &mPairSmearCiC, &b_mPairSmearCiC);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   Notify();
}

Bool_t HggOutputReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HggOutputReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HggOutputReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HggOutputReader_cxx
