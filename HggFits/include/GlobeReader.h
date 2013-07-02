//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 22 13:29:55 2013 by ROOT version 5.32/00
// from TChain spin_trees/grav2pm_m125_8TeV/
//////////////////////////////////////////////////////////

#ifndef GlobeReader_h
#define GlobeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class GlobeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           category;
   Float_t         evweight;
   Double_t        higgs_px;
   Double_t        higgs_py;
   Double_t        higgs_pz;
   Double_t        higgs_E;
   Double_t        lead_px;
   Double_t        lead_py;
   Double_t        lead_pz;
   Double_t        lead_E;
   Double_t        sublead_px;
   Double_t        sublead_py;
   Double_t        sublead_pz;
   Double_t        sublead_E;
   Double_t        costheta_cs;
   Double_t        costheta_hx;
   Double_t        lead_calo_eta;
   Double_t        lead_calo_phi;
   Float_t         lead_r9;
   Double_t        sublead_calo_eta;
   Double_t        sublead_calo_phi;
   Float_t         sublead_r9;
   Float_t         diphoton_bdt;
   Double_t        higgs_mass;
   Double_t        sigmaMrv;
   Double_t        sigmaMwv;
   Double_t        lead_sigmaE;
   Double_t        lead_sigmaE_nosmear;
   Double_t        sublead_sigmaE;
   Double_t        sublead_sigmaE_nosmear;
   Float_t         sigmaMoMrv;
   Float_t         sigmaMoMwv;
   Float_t         vtx_prob;
   Float_t         leadPtoM;
   Float_t         subleadPtoM;
   Float_t         leadEta;
   Float_t         subleadEta;
   Float_t         cosDphi;
   Float_t         lead_id_mva;
   Float_t         sublead_id_mva;

   // List of branches
   TBranch        *b_category;   //!
   TBranch        *b_evweight;   //!
   TBranch        *b_higgs_px;   //!
   TBranch        *b_higgs_py;   //!
   TBranch        *b_higgs_pz;   //!
   TBranch        *b_higgs_E;   //!
   TBranch        *b_lead_px;   //!
   TBranch        *b_lead_py;   //!
   TBranch        *b_lead_pz;   //!
   TBranch        *b_lead_E;   //!
   TBranch        *b_sublead_px;   //!
   TBranch        *b_sublead_py;   //!
   TBranch        *b_sublead_pz;   //!
   TBranch        *b_sublead_E;   //!
   TBranch        *b_costheta_cs;   //!
   TBranch        *b_costheta_hx;   //!
   TBranch        *b_lead_calo_eta;   //!
   TBranch        *b_lead_calo_phi;   //!
   TBranch        *b_lead_r9;   //!
   TBranch        *b_sublead_calo_eta;   //!
   TBranch        *b_sublead_calo_phi;   //!
   TBranch        *b_sublead_r9;   //!
   TBranch        *b_diphoton_bdt;   //!
   TBranch        *b_higgs_mass;   //!
   TBranch        *b_sigmaMrv;   //!
   TBranch        *b_sigmaMwv;   //!
   TBranch        *b_lead_sigmaE;   //!
   TBranch        *b_lead_sigmaE_nosmear;   //!
   TBranch        *b_sublead_sigmaE;   //!
   TBranch        *b_sublead_sigmaE_nosmear;   //!
   TBranch        *b_sigmaMoMrv;   //!
   TBranch        *b_sigmaMoMwv;   //!
   TBranch        *b_vtx_prob;   //!
   TBranch        *b_leadPtoM;   //!
   TBranch        *b_subleadPtoM;   //!
   TBranch        *b_leadEta;   //!
   TBranch        *b_subleadEta;   //!
   TBranch        *b_cosDphi;   //!
   TBranch        *b_lead_id_mva;   //!
   TBranch        *b_sublead_id_mva;   //!

   GlobeReader(TTree *tree=0);
   virtual ~GlobeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GlobeReader_cxx
GlobeReader::GlobeReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("spin_trees/grav2pm_m125_8TeV",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("spin_trees/grav2pm_m125_8TeV","");
      chain->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/analyzed/spinstudies/moriond2013_spin_cutbased/histograms_CMS-HGG_99.root/spin_trees/grav2pm_m125_8TeV");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

GlobeReader::~GlobeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GlobeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GlobeReader::LoadTree(Long64_t entry)
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

void GlobeReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("evweight", &evweight, &b_evweight);
   fChain->SetBranchAddress("higgs_px", &higgs_px, &b_higgs_px);
   fChain->SetBranchAddress("higgs_py", &higgs_py, &b_higgs_py);
   fChain->SetBranchAddress("higgs_pz", &higgs_pz, &b_higgs_pz);
   fChain->SetBranchAddress("higgs_E", &higgs_E, &b_higgs_E);
   fChain->SetBranchAddress("lead_px", &lead_px, &b_lead_px);
   fChain->SetBranchAddress("lead_py", &lead_py, &b_lead_py);
   fChain->SetBranchAddress("lead_pz", &lead_pz, &b_lead_pz);
   fChain->SetBranchAddress("lead_E", &lead_E, &b_lead_E);
   fChain->SetBranchAddress("sublead_px", &sublead_px, &b_sublead_px);
   fChain->SetBranchAddress("sublead_py", &sublead_py, &b_sublead_py);
   fChain->SetBranchAddress("sublead_pz", &sublead_pz, &b_sublead_pz);
   fChain->SetBranchAddress("sublead_E", &sublead_E, &b_sublead_E);
   fChain->SetBranchAddress("costheta_cs", &costheta_cs, &b_costheta_cs);
   fChain->SetBranchAddress("costheta_hx", &costheta_hx, &b_costheta_hx);
   fChain->SetBranchAddress("lead_calo_eta", &lead_calo_eta, &b_lead_calo_eta);
   fChain->SetBranchAddress("lead_calo_phi", &lead_calo_phi, &b_lead_calo_phi);
   fChain->SetBranchAddress("lead_r9", &lead_r9, &b_lead_r9);
   fChain->SetBranchAddress("sublead_calo_eta", &sublead_calo_eta, &b_sublead_calo_eta);
   fChain->SetBranchAddress("sublead_calo_phi", &sublead_calo_phi, &b_sublead_calo_phi);
   fChain->SetBranchAddress("sublead_r9", &sublead_r9, &b_sublead_r9);
   fChain->SetBranchAddress("diphoton_bdt", &diphoton_bdt, &b_diphoton_bdt);
   fChain->SetBranchAddress("higgs_mass", &higgs_mass, &b_higgs_mass);
   fChain->SetBranchAddress("sigmaMrv", &sigmaMrv, &b_sigmaMrv);
   fChain->SetBranchAddress("sigmaMwv", &sigmaMwv, &b_sigmaMwv);
   fChain->SetBranchAddress("lead_sigmaE", &lead_sigmaE, &b_lead_sigmaE);
   fChain->SetBranchAddress("lead_sigmaE_nosmear", &lead_sigmaE_nosmear, &b_lead_sigmaE_nosmear);
   fChain->SetBranchAddress("sublead_sigmaE", &sublead_sigmaE, &b_sublead_sigmaE);
   fChain->SetBranchAddress("sublead_sigmaE_nosmear", &sublead_sigmaE_nosmear, &b_sublead_sigmaE_nosmear);
   fChain->SetBranchAddress("sigmaMoMrv", &sigmaMoMrv, &b_sigmaMoMrv);
   fChain->SetBranchAddress("sigmaMoMwv", &sigmaMoMwv, &b_sigmaMoMwv);
   fChain->SetBranchAddress("vtx_prob", &vtx_prob, &b_vtx_prob);
   fChain->SetBranchAddress("leadPtoM", &leadPtoM, &b_leadPtoM);
   fChain->SetBranchAddress("subleadPtoM", &subleadPtoM, &b_subleadPtoM);
   fChain->SetBranchAddress("leadEta", &leadEta, &b_leadEta);
   fChain->SetBranchAddress("subleadEta", &subleadEta, &b_subleadEta);
   fChain->SetBranchAddress("cosDphi", &cosDphi, &b_cosDphi);
   fChain->SetBranchAddress("lead_id_mva", &lead_id_mva, &b_lead_id_mva);
   fChain->SetBranchAddress("sublead_id_mva", &sublead_id_mva, &b_sublead_id_mva);
   Notify();
}

Bool_t GlobeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GlobeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GlobeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GlobeReader_cxx
