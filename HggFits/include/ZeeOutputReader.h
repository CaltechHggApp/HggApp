//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 20 10:45:49 2012 by ROOT version 5.32/00
// from TTree ZeeOutput/
// found on file: ../../skim/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Selected.root
//////////////////////////////////////////////////////////

#ifndef ZeeOutputReader_h
#define ZeeOutputReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ZeeOutputReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         mass;
   Float_t         Ele1pt;
   Float_t         Ele1eta;
   Float_t         Ele1phi;
   Float_t         Ele1E;
   Float_t         Ele1Epho;
   Float_t         Ele2pt;
   Float_t         Ele2eta;
   Float_t         Ele2phi;
   Float_t         Ele2E;
   Float_t         Ele2Epho;
   Float_t         Ele1r9;
   Float_t         Ele1mva;
   Float_t         Ele1etaSC;
   Float_t         Ele1sigEoE;
   Float_t         Ele1sigEoEpho;
   Float_t         Ele1sigEscaleoEpho;
   Float_t         Ele2r9;
   Float_t         Ele2mva;
   Float_t         Ele2etaSC;
   Float_t         Ele2sigEoE;
   Float_t         Ele2sigEoEpho;
   Float_t         Ele2sigEscaleoEpho;
   Float_t         passloose;
   Float_t         passtight;
   Float_t         passmva;
   Float_t         nEle;

   // List of branches
   TBranch        *b_mass;   //!
   TBranch        *b_Ele1pt;   //!
   TBranch        *b_Ele1eta;   //!
   TBranch        *b_Ele1phi;   //!
   TBranch        *b_Ele1E;   //!
   TBranch        *b_Ele1Epho;   //!
   TBranch        *b_Ele2pt;   //!
   TBranch        *b_Ele2eta;   //!
   TBranch        *b_Ele2phi;   //!
   TBranch        *b_Ele2E;   //!
   TBranch        *b_Ele2Epho;   //!
   TBranch        *b_Ele1r9;   //!
   TBranch        *b_Ele1mva;   //!
   TBranch        *b_Ele1etaSC;   //!
   TBranch        *b_Ele1sigEoE;   //!
   TBranch        *b_Ele1sigEoEpho;   //!
   TBranch        *b_Ele1sigEscaleoEpho;   //!
   TBranch        *b_Ele2r9;   //!
   TBranch        *b_Ele2mva;   //!
   TBranch        *b_Ele2etaSC;   //!
   TBranch        *b_Ele2sigEoE;   //!
   TBranch        *b_Ele2sigEoEpho;   //!
   TBranch        *b_Ele2sigEscaleoEpho;   //!
   TBranch        *b_passloose;   //!
   TBranch        *b_passtight;   //!
   TBranch        *b_passmva;   //!
   TBranch        *b_nEle;   //!

   ZeeOutputReader(TTree *tree=0);
   virtual ~ZeeOutputReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ZeeOutputReader_cxx
ZeeOutputReader::ZeeOutputReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../skim/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Selected.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../skim/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Selected.root");
      }
      f->GetObject("ZeeOutput",tree);

   }
   Init(tree);
}

ZeeOutputReader::~ZeeOutputReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZeeOutputReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZeeOutputReader::LoadTree(Long64_t entry)
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

void ZeeOutputReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("Ele1pt", &Ele1pt, &b_Ele1pt);
   fChain->SetBranchAddress("Ele1eta", &Ele1eta, &b_Ele1eta);
   fChain->SetBranchAddress("Ele1phi", &Ele1phi, &b_Ele1phi);
   fChain->SetBranchAddress("Ele1E", &Ele1E, &b_Ele1E);
   fChain->SetBranchAddress("Ele1Epho", &Ele1Epho, &b_Ele1Epho);
   fChain->SetBranchAddress("Ele2pt", &Ele2pt, &b_Ele2pt);
   fChain->SetBranchAddress("Ele2eta", &Ele2eta, &b_Ele2eta);
   fChain->SetBranchAddress("Ele2phi", &Ele2phi, &b_Ele2phi);
   fChain->SetBranchAddress("Ele2E", &Ele2E, &b_Ele2E);
   fChain->SetBranchAddress("Ele2Epho", &Ele2Epho, &b_Ele2Epho);
   fChain->SetBranchAddress("Ele1r9", &Ele1r9, &b_Ele1r9);
   fChain->SetBranchAddress("Ele1mva", &Ele1mva, &b_Ele1mva);
   fChain->SetBranchAddress("Ele1etaSC", &Ele1etaSC, &b_Ele1etaSC);
   fChain->SetBranchAddress("Ele1sigEoE", &Ele1sigEoE, &b_Ele1sigEoE);
   fChain->SetBranchAddress("Ele1sigEoEpho", &Ele1sigEoEpho, &b_Ele1sigEoEpho);
   fChain->SetBranchAddress("Ele1sigEscaleoEpho", &Ele1sigEscaleoEpho, &b_Ele1sigEscaleoEpho);
   fChain->SetBranchAddress("Ele2r9", &Ele2r9, &b_Ele2r9);
   fChain->SetBranchAddress("Ele2mva", &Ele2mva, &b_Ele2mva);
   fChain->SetBranchAddress("Ele2etaSC", &Ele2etaSC, &b_Ele2etaSC);
   fChain->SetBranchAddress("Ele2sigEoE", &Ele2sigEoE, &b_Ele2sigEoE);
   fChain->SetBranchAddress("Ele2sigEoEpho", &Ele2sigEoEpho, &b_Ele2sigEoEpho);
   fChain->SetBranchAddress("Ele2sigEscaleoEpho", &Ele2sigEscaleoEpho, &b_Ele2sigEscaleoEpho);
   fChain->SetBranchAddress("passloose", &passloose, &b_passloose);
   fChain->SetBranchAddress("passtight", &passtight, &b_passtight);
   fChain->SetBranchAddress("passmva", &passmva, &b_passmva);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   Notify();
}

Bool_t ZeeOutputReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZeeOutputReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZeeOutputReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZeeOutputReader_cxx
