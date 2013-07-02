#include <ZeeSelector.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"
using namespace std;
#define PI atan(1.0)*4
ZeeSelector::ZeeSelector():
  fChain(0),
  valid(false),
  isData_(true),
  Electrons(0)
{
}

ZeeSelector::ZeeSelector(vector<string> fNames, string treeName,string outFName):
  isData_(true),
  Electrons(0)
{
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

ZeeSelector::~ZeeSelector(){
  delete fChain;
}


void ZeeSelector::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int ZeeSelector::init(){
  if(!valid) return -1;

  elecorr = new HggEGEnergyCorrector(0,config,isData_);
  
  this->setBranchAddresses();
  this->setupOutputTree();

  ReadConfig cfg;
  int retcode = cfg.read(config);
  if(retcode != 0){
    cout << "Error reading configuration file!" <<std::endl
	 << "Error Code: " << retcode << std::endl;
    valid = false;
    return -1;
  }

  doScale      = cfg.getParameter("redoScale").compare("yes")==0;
  doSmear      = cfg.getParameter("redoSmear").compare("yes")==0;
  
  if(doScale){
    std::string scaleCFG = cfg.getParameter("EnergyScaleCFG");
    scale = new HggEnergyScale(scaleCFG);
  }
  if(doSmear){
    std::string smearCFG = cfg.getParameter("EnergySmearCFG");
    smear = new HggEnergyScale(smearCFG);    
    applyScaleSmear = atoi(cfg.getParameter("ScaleSmear").c_str());
    std::cout << "Doing smearing.  Config: " << smearCFG << "  ScaleSmear: " << applyScaleSmear <<std::endl;
  }

}

void swap(VecbosEle& e1, VecbosEle &e2){
  VecbosEle t= e1; e1=e2; e2=e1;
}

void ZeeSelector::Loop(){
  if(!valid) return;

  this->init();
  cout << "Getting Entries ... "  << endl;
  int nEntries = fChain->GetEntries();
  int nbytes = 0;

  for (int eventi=0; eventi<nEntries; eventi++) {
    if(eventi%500==0) {
      cout << ">> Processing Entry " << eventi << "/" << nEntries << endl;
    }    
    nbytes += fChain->GetEvent(eventi);

    // Reset mass reference
    DZmassref = 100;
    isZmass = false;
    // Require the event to have at least two electrons
    if (nEle > 1) {   

      //apply scale and smear to the electrons
      for(int i=0;i<nEle;i++){
	VecbosEle ele = Electrons->at(i);
	if(doScale){
	  std::pair<float,float> dE = scale->getDEoE(ele,runNumber);
	  ele.dEoE    = dE.first;
	  ele.dEoEErr = 0;
	  ele.scaledEnergy = ele.correctedEnergy*(ele.dEoE);
	  ele.scaledEnergyError = ele.correctedEnergyError*((ele.dEoE+ele.dEoEErr));
	}
	if(doSmear){
	  ele.scaledEnergy = ele.correctedEnergy;
	  ele.scaledEnergyError = ele.correctedEnergyError;
	  std::pair<float,float> dE = smear->getDEoE(*pho,applyScaleSmear);
	  ele.dEoE    = dE.first;
	  ele.dEoEErr = dE.second;
	}
      }


      // Choose First Electron
      for (int k=0; k<nEle-1;k++) { 
	VecbosEle ele1 = Electrons->at(k);
	if(ele1.SC.index<0) continue;
	if(!passPresel(ele1)) continue;
	// Choose Second Electron
	for (int j=k+1; j<nEle;j++) {   
	  lpass   = 0;
	  tpass   = 0;
	  mvapass = 0;
	  Zeemass = 0;
	  VecbosEle ele2 = Electrons->at(j);
	  if(ele2.SC.index<0) continue;
	  if(!passPresel(ele2)) continue;
	  // Verify Opposite Charges
	  if (ele1.charge != ele2.charge) {

	    if(ele1.pt < ele2.pt) swap(ele1,ele2);

	    // Calculate Invariant Mass for Z Boson from Two Electrons
	    TLorentzVector Ele1;
	    Ele1.SetPtEtaPhiM(ele1.correctedEnergy/cosh(ele1.eta),ele1.eta,ele1.phi,0);
	    TLorentzVector Ele2;
	    Ele2.SetPtEtaPhiM(ele2.correctedEnergy/cosh(ele2.eta),ele2.eta,ele2.phi,0);
	    TLorentzVector Zee = Ele1 + Ele2;
	    Zeemass = Zee.M();

	    //preselection cuts:
	    

	    PFIsoOverPT1 = (ele1.dr03ChargedHadronPFIso + max(0.d,(ele1.dr03NeutralHadronPFIso+ele1.dr03PhotonPFIso) - pow(0.3, 2)*PI*max(0.f,rho)))/(ele1.pt);
	    PFIsoOverPT2 = (ele2.dr03ChargedHadronPFIso + max(0.d,(ele2.dr03NeutralHadronPFIso+ele2.dr03PhotonPFIso) - pow(0.3, 2)*PI*max(0.f,rho)))/(ele2.pt);
	    
	    // Loose Cuts - WP 90
	    if (ele1.correctedEnergy>25 && ele2.correctedEnergy>25) {
	      if (((fabs(ele1.SC.eta)<1.44 && fabs(ele1.dEtaSCTrackAtVtx)<0.007 && fabs(ele1.dPhiSCTrackAtVtx)<0.15 && ele1.SC.sigmaIEtaIEta<0.01 && ele1.SC.HoverE<0.12) || (fabs(ele1.SC.eta)>1.52 && fabs(ele1.dEtaSCTrackAtVtx)<0.009 && fabs(ele1.dPhiSCTrackAtVtx)<0.1 && ele1.SC.sigmaIEtaIEta<0.03 && ele1.SC.HoverE<0.10)) && ((fabs(ele2.SC.eta)<1.44 && fabs(ele2.dEtaSCTrackAtVtx)<0.007 && fabs(ele2.dPhiSCTrackAtVtx)<0.15 && ele2.SC.sigmaIEtaIEta<0.01 && ele2.SC.HoverE<0.12) || (fabs(ele2.SC.eta)>1.52 && fabs(ele2.dEtaSCTrackAtVtx)<0.009 && fabs(ele2.dPhiSCTrackAtVtx)<0.1 && ele2.SC.sigmaIEtaIEta<0.03 && ele2.SC.HoverE<0.10))) {
		if (PFIsoOverPT1<0.15 && PFIsoOverPT2<0.15 && ele1.hasMatchedConversion == false && ele2.hasMatchedConversion == false && ele1.expInnerLayersHits!=-999 && ele2.expInnerLayersHits!=-999) {
		  if (ele1.d0Track != 999 && ele2.d0Track != 999 && ele1.dzTrack != 999 && ele2.dzTrack != 999) {
		    lpass = 1.0;
		    
		    // Tight Cuts - WP 70
		    if (((fabs(ele1.SC.eta)<1.44 && fabs(ele1.dEtaSCTrackAtVtx)<.004 && fabs(ele1.dPhiSCTrackAtVtx)<0.03) || (fabs(ele1.SC.eta)>1.52 && fabs(ele1.dEtaSCTrackAtVtx)<0.005 && fabs(ele1.dPhiSCTrackAtVtx)<0.02)) && ((fabs(ele2.SC.eta)<1.44 && fabs(ele1.dEtaSCTrackAtVtx)<0.004 && fabs(ele2.dPhiSCTrackAtVtx)<0.03) || (fabs(ele2.SC.eta)>1.52 && fabs(ele2.dEtaSCTrackAtVtx)<0.005 && fabs(ele2.dPhiSCTrackAtVtx)<0.02))) {
		      if (PFIsoOverPT1<0.10 && PFIsoOverPT2<0.10) {
			tpass = 1.0;
		      };
		    };
		  };
		};
	      };
	    
		    
	      // MVA ID Cuts
	      if (((fabs(ele1.SC.eta)<0.8 && ele1.idMVA>0.5) || (fabs(ele1.SC.eta)>0.8 && fabs(ele1.SC.eta)<1.479 && ele1.idMVA>0.120) || (fabs(ele1.SC.eta)>1.479 && ele1.idMVA>0.6)) && ((fabs(ele2.SC.eta)<0.8 && ele2.idMVA>0.5) || (fabs(ele2.SC.eta)>0.8 && fabs(ele2.SC.eta)<1.479 && ele2.idMVA>0.120) || (fabs(ele2.SC.eta)>1.479 && ele2.idMVA>0.6))) {
		mvapass = 1.0;
	      };
	    };
	  	    
	    // Calculate Difference from True Z Mass 
	    DZmass = fabs(Zeemass - 91.2);
	    // Compare the proximity of uncut Z mass to real Z mass with other electron pairs in event
	    if (DZmass < DZmassref && (lpass + mvapass + tpass)>=0) {
	      elecorr->setEventInfo(rho,nVtx);

	      std::pair<double,double> ele1corr = elecorr->electronEnergyCorrector_May2012(ele1,false);
	      std::pair<double,double> ele2corr = elecorr->electronEnergyCorrector_May2012(ele2,false);
	      std::pair<double,double> ele1corrScale = elecorr->electronEnergyCorrector_May2012(ele1,true);
	      std::pair<double,double> ele2corrScale = elecorr->electronEnergyCorrector_May2012(ele2,true);

		// Reset the selected Z mass and reference point to this pair
		mass = Zeemass;
		DZmassref = DZmass;          

		Ele1pt = ele1.pt;
		Ele1eta = ele1.eta;
		Ele1phi = ele1.phi;
		Ele1E   = ele1.correctedEnergy;
		Ele1Epho = ele1corr.first;
		Ele1sigEoE = ele1.correctedEnergyError/ele1.correctedEnergy;
		Ele1sigEoEpho = ele1corr.second/ele1corr.first;
		Ele1sigEscaleoEpho = ele1corrScale.second/ele1corrScale.first;
		
		Ele2pt = ele2.pt;
		Ele2eta = ele2.eta;
		Ele2phi = ele2.phi;
		Ele2E   = ele2.correctedEnergy;
		Ele2Epho = ele2corr.first;
		Ele2sigEoE = ele2.correctedEnergyError/ele2.correctedEnergy;
		Ele2sigEoEpho = ele2corr.second/ele2corr.first;
		Ele2sigEscaleoEpho = ele2corrScale.second/ele2corrScale.first;

		Ele1etaSC = ele1.SC.eta;
		Ele2etaSC = ele2.SC.eta;
		Ele1r9  = ele1.SC.r9;
		Ele2r9  = ele2.SC.r9;
		passloose = lpass;
		passtight = tpass;
		passmva   = mvapass;
		Ele1mva = ele1.idMVA;
		Ele2mva = ele2.idMVA;
		isZmass = true;
	    };
	  };
	};
      };
    };
    if (isZmass ==  true) {
      nEleOut = nEle;
      outTree->Fill();
    };
  };

  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outTree->Write();
  f->Close();
  std::map<std::string,TH1F*>::iterator it;
}



void ZeeSelector::setBranchAddresses(){
  if(!valid) return;
  //Event info
  fChain->SetBranchAddress("runNumber",&runNumber);
  fChain->SetBranchAddress("evtNumber",&evtNumber);
  //fChain->SetBranchAddress("isRealData",&_isData);
  
  //information for the Electrons
  fChain->SetBranchAddress("nEle",&nEle);
  fChain->SetBranchAddress("Electrons",&Electrons);
 
  // information for the vertex
  fChain->SetBranchAddress("nVtx", &nVtx);
  fChain->SetBranchAddress("rho", &rho);

}

void ZeeSelector::setupOutputTree(){
  outTree = new TTree("ZeeOutput","");
  outTree->Branch("mass",&mass,"mass");

  outTree->Branch("Ele1pt",&Ele1pt,"Ele1pt");
  outTree->Branch("Ele1eta",&Ele1eta,"Ele1eta");
  outTree->Branch("Ele1phi",&Ele1phi,"Ele1phi");
  outTree->Branch("Ele1E",&Ele1E,"Ele1E");
  outTree->Branch("Ele1Epho",&Ele1Epho,"Ele1Epho");

  outTree->Branch("Ele2pt",&Ele2pt,"Ele2pt");
  outTree->Branch("Ele2eta",&Ele2eta,"Ele2eta");
  outTree->Branch("Ele2phi",&Ele2phi,"Ele2phi");
  outTree->Branch("Ele2E",&Ele2E,"Ele2E");
  outTree->Branch("Ele2Epho",&Ele2Epho,"Ele2Epho");

  outTree->Branch("Ele1r9",&Ele1r9,"Ele1r9");
  outTree->Branch("Ele1mva",&Ele1mva,"Ele1mva");
  outTree->Branch("Ele1etaSC",&Ele1etaSC,"Ele1etaSC");
  outTree->Branch("Ele1sigEoE",&Ele1sigEoE,"Ele1sigEoE");
  outTree->Branch("Ele1sigEoEpho",&Ele1sigEoEpho,"Ele1sigEoEpho");
  outTree->Branch("Ele1sigEscaleoEpho",&Ele1sigEscaleoEpho,"Ele1sigEscaleoEpho");
  outTree->Branch("Ele2r9",&Ele2r9,"Ele2r9");
  outTree->Branch("Ele2mva",&Ele2mva,"Ele2mva");
  outTree->Branch("Ele2etaSC",&Ele2etaSC,"Ele2etaSC");
  outTree->Branch("Ele2sigEoE",&Ele2sigEoE,"Ele2sigEoE");
  outTree->Branch("Ele2sigEoEpho",&Ele2sigEoEpho,"Ele2sigEoEpho");
  outTree->Branch("Ele2sigEscaleoEpho",&Ele2sigEscaleoEpho,"Ele2sigEscaleoEpho");

  outTree->Branch("passloose",&passloose,"passloose");  
  outTree->Branch("passtight",&passtight,"passtight");
  outTree->Branch("passmva",&passmva,"passmva");
  outTree->Branch("nEle",&nEleOut,"nEle");
}

bool ZeeSelector::passPresel(VecbosEle &ele){
  //cuts barrel_low,endcap_low,barrel_high,endcap_high
  float ecal[4] = {4.,4.,50.,50.};
  float hcal[4] = {4.,4.,50.,50.};
  float trk[4] = {4.,4.,50.,50.};
  float he[4] = {0.075,0.075,0.082,0.075};
  float sieie[4] = { 0.014,0.034,0.014,0.034};

  int index = (ele.SC.r9>0.9)*2+(fabs(ele.SC.eta)>1.48);

  if( ele.dr03EcalRecHitSumEt - 0.012*ele.pt > ecal[index] ||
      ele.dr03HcalTowerSumEt  - 0.005*ele.pt > hcal[index] ||
      ele.dr03TkSumPt         - 0.002*ele.pt > trk[index]) return false;
  if( ele.HoverE > he[index] || ele.SC.sigmaIEtaIEta > sieie[index]) return false;
  
  return true;
  
}






