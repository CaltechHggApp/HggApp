#include <VecbosEGObject.hh>
#include <src/HggPhysUtils.cc>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <iostream>
#include <TH2D.h>

int main(int argc, char* argv[]){
  char inputFileName[150];
  if( argc < 2 ) return 0;
  strcpy(inputFileName,argv[1]);

  TChain *tree= new TChain("HggReduce");
  tree->AddFile(inputFileName);

  std::vector<VecbosPho> *Photons;
  Int_t nPho;
  const int MAXGenSaved = 1000;
  //gen-leve phton 
  int nGenPho;
  float etaGenPho[MAXGenSaved];
  float phiGenPho[MAXGenSaved];
  float ptGenPho[MAXGenSaved];
  float energyGenPho[MAXGenSaved];
  int pidMomGenPho[MAXGenSaved];
  int indMomGenPho[MAXGenSaved];
  int statusGenPho[MAXGenSaved];
  int indexGenPho[MAXGenSaved];

  tree->SetBranchAddress("Photons",&Photons);
  tree->SetBranchAddress("nPho",&nPho); 
  
  tree->SetBranchAddress("nGenPho",&nGenPho);
  tree->SetBranchAddress("etaGenPho",etaGenPho);
  tree->SetBranchAddress("phiGenPho",phiGenPho);
  tree->SetBranchAddress("ptGenPho",ptGenPho);
  tree->SetBranchAddress("energyGenPho",energyGenPho);
  
  TH2D* ErecoCorEtrue = new TH2D("ERecoCorByETrue","",100,-5,5,100,0,2);
  TH2D* ErecoEtrue = new TH2D("ERecoByETrue","",100,-5,5,100,0,2);
  TH1D* EBHighR9   = new TH1D("EBHighR9","",100,0,2);
  TH1D* EBLowR9   = new TH1D("EBLowR9","",100,0,2);
  TH1D* EEHighR9   = new TH1D("EEHighR9","",100,0,2);
  TH1D* EELowR9   = new TH1D("EELowR9","",100,0,2);

  TH1D* EBHighCorR9   = new TH1D("EBHighCorR9","",100,0,2);
  TH1D* EBLowCorR9   = new TH1D("EBLowCorR9","",100,0,2);
  TH1D* EEHighCorR9   = new TH1D("EEHighCorR9","",100,0,2);
  TH1D* EELowCorR9   = new TH1D("EELowCorR9","",100,0,2);

  EBHighR9->SetXTitle("E_{reco}/E_{true}");
  EBLowR9->SetXTitle("E_{reco}/E_{true}");
  EEHighR9->SetXTitle("E_{reco}/E_{true}");
  EELowR9->SetXTitle("E_{reco}/E_{true}");

  EBHighCorR9->SetXTitle("E_{reco}/E_{true}");
  EBLowCorR9->SetXTitle("E_{reco}/E_{true}");
  EEHighCorR9->SetXTitle("E_{reco}/E_{true}");
  EELowCorR9->SetXTitle("E_{reco}/E_{true}");

  const float maxDR=0.2;
  Int_t iEntry=-1;
  while(tree->GetEntry(++iEntry)){
    if(iEntry%100==0) std::cout << "Processing Entry: " << iEntry << std::endl;
    int j=0;
    std::vector<VecbosPho>::const_iterator recoI;
    for(j=0,recoI = Photons->begin(); j<nPho;recoI++,j++){
      if(recoI->energy/TMath::CosH(recoI->eta) < 25) continue;
      double minDR=9999;
      int minIndex=-1;
      for(int i=0;i<nGenPho;i++){
	double dr = DeltaR(recoI->SC.eta,etaGenPho[i],recoI->SC.phi,phiGenPho[i]);
	if(dr<minDR && dr<maxDR){
	  minDR = dr;
	  minIndex=i;
	}
      }
      if(minIndex==-1) continue;
      ErecoEtrue->Fill(recoI->SC.eta,recoI->energy/energyGenPho[minIndex]);
      ErecoCorEtrue->Fill(recoI->SC.eta,recoI->finalEnergy/energyGenPho[minIndex]);      
      double r9 = recoI->SC.e3x3/recoI->SC.rawE*recoI->SC.r9Scale;
      if(fabs(recoI->SC.eta) < 1.48 && r9 >= 0.94) EBHighR9->Fill(recoI->energy/energyGenPho[minIndex]);
      if(fabs(recoI->SC.eta) < 1.48 && r9 < 0.94) EBLowR9->Fill(recoI->energy/energyGenPho[minIndex]);
      if(fabs(recoI->SC.eta) >= 1.48 && r9 >= 0.94) EEHighR9->Fill(recoI->energy/energyGenPho[minIndex]);
      if(fabs(recoI->SC.eta) >= 1.48 && r9 < 0.94) EELowR9->Fill(recoI->energy/energyGenPho[minIndex]);

      if(fabs(recoI->SC.eta) < 1.48 && r9 >= 0.94) EBHighCorR9->Fill(recoI->finalEnergy/energyGenPho[minIndex]);
      if(fabs(recoI->SC.eta) < 1.48 && r9 < 0.94) EBLowCorR9->Fill(recoI->finalEnergy/energyGenPho[minIndex]);
      if(fabs(recoI->SC.eta) >= 1.48 && r9 >= 0.94) EEHighCorR9->Fill(recoI->finalEnergy/energyGenPho[minIndex]);
      if(fabs(recoI->SC.eta) >= 1.48 && r9 < 0.94) EELowCorR9->Fill(recoI->finalEnergy/energyGenPho[minIndex]);

    }//for(recoI = Photons->begin()...
  }//while(tree->GetEntry(++iEntry))
  std::cout << "Writing..." << std::endl;
  TFile *out = TFile::Open("ERecoOverETrue.root","RECREATE");
  out->cd();
  ErecoEtrue->Write();
  ErecoCorEtrue->Write();
  EBHighR9->Write();
  EBLowR9->Write();
  EEHighR9->Write();
  EELowR9->Write();

  EBHighCorR9->Write();
  EBLowCorR9->Write();
  EEHighCorR9->Write();
  EELowCorR9->Write();

  out->Close();
  return 0;
}
