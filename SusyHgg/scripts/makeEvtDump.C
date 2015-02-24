#include "TTree.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "TTreeFormula.h"
#include "assert.h"

void makeEvtDump(TTree* SusyHggTree,const char* fileName) {

  int run,lumi,evt,Njets;
  float MR, Rsq,ucRsq, nj, pho1_pt,pho1_eta,pho1_phi,pho2_pt,pho2_eta,pho2_phi,mgg,ptgg;

  SusyHggTree->SetBranchAddress("run",&run); 
 SusyHggTree->SetBranchAddress("lumi",&lumi);
  SusyHggTree->SetBranchAddress("evt",&evt);
 SusyHggTree->SetBranchAddress("MR",&MR);
 SusyHggTree->SetBranchAddress("Rsq",&ucRsq);
 SusyHggTree->SetBranchAddress("t1Rsq",&Rsq);

 SusyHggTree->SetBranchAddress("Njets",&Njets);

  SusyHggTree->SetBranchAddress("pho1_pt",&pho1_pt); SusyHggTree->SetBranchAddress("pho1_eta",&pho1_eta); SusyHggTree->SetBranchAddress("pho1_phi",&pho1_phi);
  SusyHggTree->SetBranchAddress("pho2_pt",&pho2_pt); SusyHggTree->SetBranchAddress("pho2_eta",&pho2_eta); SusyHggTree->SetBranchAddress("pho2_phi",&pho2_phi);
  SusyHggTree->SetBranchAddress("mgg",&mgg); SusyHggTree->SetBranchAddress("ptgg",&ptgg);

  ofstream output(fileName);
  output << "run/I:evt:lumi:MR/F:Rsq:uncorrRsq:mgg/F:ptgg:pho1_pt:pho1_eta:pho1_phi:pho2_pt:pho2_eta:pho2_phi\n";
  Long64_t iEntry=-1;
  while(SusyHggTree->GetEntry(++iEntry)){ 
    if(evt==480229083) std::cout << run << " " << evt << std::endl;
    output << run << " " << evt << " " << lumi << " " << MR << " " << Rsq << " " << ucRsq << " " << mgg << " " << ptgg << " " << pho1_pt << " " << pho1_eta << " " << pho1_phi << " " << pho2_pt << " " << pho2_eta << " " << pho2_phi << std::endl; 
  }

  output.close();
}

TTree* createCopy(TTree* orig, TTree** friends, TString* selections,int nFriends) {
  Long64_t iEntry=-1;

  TTree* copy = orig->CloneTree(0);
  std::vector<TTreeFormula*> forms;
  for(int i=0;i<nFriends;i++) {
    forms.push_back(new TTreeFormula(Form("form_%d",i),selections[i],friends[i]));
  }
  int nPass=0;
  while(orig->GetEntry(++iEntry)) {
    bool pass = true;
    for(int i=0;i<nFriends;i++) {
      friends[i]->GetEntry(iEntry);
      pass = pass && (forms.at(i)->EvalInstance()>0);
    }    
    if(pass) {
      copy->Fill();
      nPass++;
    }
  }
  for(int i=0;i<nFriends;i++) delete forms.at(i);
  std::cout << "nPass: " << nPass << std::endl;
  return copy;
}

void dumpAllBoxes() {
  TFile *_file0 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7.root");

  TTree* all = ((TTree*)_file0->Get("SusyHggTree"));
  //all->AddFriend("SusyHggTriggerTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__HggTRIGGER.root");
  
  TFile *_file1 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7__TRIGGER.root");
  TTree* trigger = (TTree*)_file1->Get("SusyHggTriggerTree");
  TFile *_file2 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7__Noise.root");
  TTree* noise = (TTree*)_file2->Get("SusyHggNoiseTree");
  //all->AddFriend("SusyHggTriggerTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__TRIGGER.root");
  //all->AddFriend("SusyHggNoiseTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__Noise.root");
  TTree* friends[] = {trigger,noise};
  TString selections[] = {"passTrigger","passNoise"};
  TFile *f = new TFile("/wntmp/scratch/tmp_evtdump.root","RECREATE");
  TTree* passAll = createCopy(all,friends,selections,2);

  assert(passAll->GetEntries() > 1);

  //all->Scan("run:evt:passNoise:passTrigger:pho1_pt:pho2_pt:pho1_pass_iso:pho2_pass_iso:pho1_eta:pho2_eta","evt==480229083");
  passAll->Scan("run:evt","evt==347257666 || evt==50763554 || evt==76154706");
  
  TTree *sel = passAll->CopyTree("pho1_pt>40 && pho2_pt>25 && !(pho1_pass_iso && pho2_pass_iso) && (pho1_pass_iso || pho2_pass_iso) && abs(pho1_eta)<1.5 && abs(pho2_eta)<1.5");

  TFile *f1 = new TFile("/wntmp/scratch/tmp_evtdump_HighPt.root","RECREATE");
  TTree *HighPt = sel->CopyTree("ptgg>=110 ");
  makeEvtDump(HighPt,"HighPt_evtList_RunABCD_HLT3_MetFilter2.txt");
  delete HighPt;
  f1->Close();

  TFile *f2 = new TFile("/wntmp/scratch/tmp_evtdump_Hbb.root","RECREATE");
  TTree *Hbb = sel->CopyTree("ptgg < 110 && (abs(mbb_NearH-125)<15) ");
  makeEvtDump(Hbb,"Hbb_evtList_RunABCD_HLT3_MetFilter2.txt");
  delete Hbb;
  f2->Close();

  TFile *f3 = new TFile("/wntmp/scratch/tmp_evtdump_Zbb.root","RECREATE");
  TTree *Zbb = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && (abs(mbb_NearZ-91.2)<15) ");
  makeEvtDump(Zbb,"Zbb_evtList_RunABCD_HLT3_MetFilter2.txt");
  delete Zbb;
  f3->Close();

  TFile *f4 = new TFile("/wntmp/scratch/tmp_evtdump_HighRes.root","RECREATE");
  TTree *HighRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && (pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  makeEvtDump(HighRes,"HighRes_evtList_RunABCD_HLT3_MetFilter2.txt");
  delete HighRes;
  f4->Close();

  TFile *f5 = new TFile("/wntmp/scratch/tmp_evtdump_LowRes.root","RECREATE");
  TTree *LowRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && !(pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  makeEvtDump(LowRes,"LowRes_evtList_RunABCD_HLT3_MetFilter2.txt");
  delete LowRes;
  f5->Close();

  delete sel;
  f->Close();
}


void splitAllBoxes(bool invertIso=false) {
  TFile *_file0 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7.root");

  TTree* all = ((TTree*)_file0->Get("SusyHggTree"));
  //all->AddFriend("SusyHggTriggerTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__HggTRIGGER.root");
  
  TFile *_file1 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7__TRIGGER.root");
  TTree* trigger = (TTree*)_file1->Get("SusyHggTriggerTree");
  TFile *_file2 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7__Noise.root");
  TTree* noise = (TTree*)_file2->Get("SusyHggNoiseTree");
  //all->AddFriend("SusyHggTriggerTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__TRIGGER.root");
  //all->AddFriend("SusyHggNoiseTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__Noise.root");
  TTree* friends[] = {trigger,noise};
  TString selections[] = {"passTrigger","passNoise"};
  TFile *f = new TFile("/wntmp/scratch/tmp_evtdump.root","RECREATE");
  TTree* passAll = createCopy(all,friends,selections,2);

  assert(passAll->GetEntries() > 1);
  
  TFile *ff = new TFile("/wntmp/scratch/tmp_evtdump.root","RECREATE");
  TTree *sel = 0;
  if(invertIso) {
    sel = passAll->CopyTree("pho1_pt>40 && pho2_pt>25 && !(pho1_pass_iso && pho2_pass_iso) && (pho1_pass_iso || pho2_pass_iso) && abs(pho1_eta)<1.5 && abs(pho2_eta)<1.5");
  } else {
    sel = passAll->CopyTree("pho1_pt>40 && pho2_pt>25 && pho1_pass_iso && pho2_pass_iso && abs(pho1_eta)<1.5 && abs(pho2_eta)<1.5");
  }

  TString tag = (invertIso? "_invertIso":"");
  TFile *f1 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v7"+tag+"__HighPt.root","RECREATE");
  TTree *HighPt = sel->CopyTree("ptgg>=110");
  HighPt->Write();
  f1->Close();

  TFile *f2 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v7"+tag+"__Hbb.root","RECREATE");
  TTree *Hbb = sel->CopyTree("ptgg < 110 && (abs(mbb_NearH-125)<15)");
  Hbb->Write();
  f2->Close();

  TFile *f3 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v7"+tag+"__Zbb.root","RECREATE");
  TTree *Zbb = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && (abs(mbb_NearZ-91.2)<15)");
  Zbb->Write();
  f3->Close();

  TFile *f4 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v7"+tag+"__HighRes.root","RECREATE");
  TTree *HighRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && (pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  HighRes->Write();
  f4->Close();

  TFile *f5 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v7"+tag+"__LowRes.root","RECREATE");
  TTree *LowRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && !(pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  LowRes->Write();
  f5->Close();

  delete sel;
  f->Close();
}


void splitAllBoxes_v8() {
  TFile *_file0 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v8.root");

  TTree* all = ((TTree*)_file0->Get("SusyHggTree"));
  //all->AddFriend("SusyHggTriggerTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__HggTRIGGER.root");
  
  //TFile *_file1 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7__TRIGGER.root");
  //TTree* trigger = (TTree*)_file1->Get("SusyHggTriggerTree");
  //TFile *_file2 = TFile::Open("/home/amott/raid4/DoublePhoton_22Jan2013_Run2012ABCD_v7__Noise.root");
  //TTree* noise = (TTree*)_file2->Get("SusyHggNoiseTree");
  //all->AddFriend("SusyHggTriggerTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__TRIGGER.root");
  //all->AddFriend("SusyHggNoiseTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__Noise.root");
  //TTree* friends[] = {trigger,noise};
  //TString selections[] = {"passTrigger","passNoise"};
  TFile *f = new TFile("/wntmp/scratch/tmp_evtdump.root","RECREATE");
  //TTree* passAll = createCopy(all,friends,selections,2);

  TTree* passAll=all;

  assert(passAll->GetEntries() > 1);
  
  TFile *ff = new TFile("/wntmp/scratch/tmp_evtdump.root","RECREATE");
  TTree *sel = passAll->CopyTree("ptgg>20 && pho1_pt>40 && pho2_pt>25 && pho1_pass_iso && pho2_pass_iso && abs(pho1_eta)<1.5 && abs(pho2_eta)<1.5 && mgg>103 && mgg < 500");

  TFile *f1 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v8__HighPt.root","RECREATE");
  TTree *HighPt = sel->CopyTree("ptgg>=110");
  HighPt->Write();
  f1->Close();

  TFile *f2 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v8__Hbb.root","RECREATE");
  TTree *Hbb = sel->CopyTree("ptgg < 110 && (abs(mbb_NearH-125)<15)");
  Hbb->Write();
  f2->Close();

  TFile *f3 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v8__Zbb.root","RECREATE");
  TTree *Zbb = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && (abs(mbb_NearZ-91.2)<15)");
  Zbb->Write();
  f3->Close();

  TFile *f4 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v8__HighRes.root","RECREATE");
  TTree *HighRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && (pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  HighRes->Write();
  f4->Close();

  TFile *f5 = new TFile("./DoublePhoton_22Jan2013_Run2012ABCD_v8__LowRes.root","RECREATE");
  TTree *LowRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && !(pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  LowRes->Write();
  f5->Close();

  delete sel;
  f->Close();
}


void splitAllBoxesGeneric(TString fileName, TString outputTag) {
  TFile *_file0 = TFile::Open(fileName);

  TTree* all = ((TTree*)_file0->Get("SusyHggTree"));
  //all->AddFriend("SusyHggTriggerTree","/home/amott/raid4/DoublePhoton_22Jan2013_Run2012A_v7__HggTRIGGER.root");
  
  TFile *f = new TFile("/wntmp/scratch/tmp_evtdump.root","RECREATE");
  TTree* passAll = all;

  assert(passAll->GetEntries() > 1);
  
  TFile *ff = new TFile("/wntmp/scratch/tmp_evtdump.root","RECREATE");
  TTree *sel = passAll->CopyTree("pho1_pt>40 && pho2_pt>25 && pho1_pass_iso && pho2_pass_iso && abs(pho1_eta)<1.5 && abs(pho2_eta)<1.5 && mgg>103 && mgg < 180");

  TFile *f1 = new TFile("./"+outputTag+"__HighPt.root","RECREATE");
  TTree *HighPt = sel->CopyTree("ptgg>=110");
  HighPt->Write();
  f1->Close();

  TFile *f2 = new TFile("./"+outputTag+"__Hbb.root","RECREATE");
  TTree *Hbb = sel->CopyTree("ptgg < 110 && (abs(mbb_NearH-125)<15)");
  Hbb->Write();
  f2->Close();

  TFile *f3 = new TFile("./"+outputTag+"__Zbb.root","RECREATE");
  TTree *Zbb = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && (abs(mbb_NearZ-91.2)<15)");
  Zbb->Write();
  f3->Close();

  TFile *f4 = new TFile("./"+outputTag+"__HighRes.root","RECREATE");
  TTree *HighRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && (pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  HighRes->Write();
  f4->Close();

  TFile *f5 = new TFile("./"+outputTag+"__LowRes.root","RECREATE");
  TTree *LowRes = sel->CopyTree("ptgg < 110 && !(abs(mbb_NearH-125)<15) && !(abs(mbb_NearZ-91.2)<15) && !(pho1_sigEoE<0.015 && pho2_sigEoE<0.015)");
  LowRes->Write();
  f5->Close();

  delete sel;
  f->Close();
}
