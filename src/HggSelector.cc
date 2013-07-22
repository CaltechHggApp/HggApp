#include <HggSelector.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"

#include "assert.h"
using namespace std;
using namespace TMVA;

#define debugSelector 0


void HggSelector::firstInit(){
  setProcessCollection("Muons",false);
  setProcessCollection("Electrons",false);
}

//
// Process the configuration file(s)
//
void HggSelector::processConfig(ReadConfig &cfg){
  weightFile_diPho = cfg.getParameter("weightFile_diPho");
  methodName_diPho  = cfg.getParameter("methodName_diPho");
  triggers = cfg.getTokens("Triggers",",");
  doRegression = cfg.getParameter("redoRegression").compare("yes")==0;
  doScale      = cfg.getParameter("redoScale").compare("yes")==0;
  doSmear      = cfg.getParameter("redoSmear").compare("yes")==0;

  if(doRegression){
    //corrector = new HggEGEnergyCorrector(this,
  }
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

  cout << "Parameters: " << endl
       << weightFile_diPho << endl
       << methodName_diPho << endl;

  string MassResConfig = cfg.getParameter("MassResolutionConfig");
  massRes = new HggMassResolution();
  massRes->setConfigFile(MassResConfig);
  massRes->init();

  //get kinematic cuts
  string leadPhoEtMinS  =    cfg.getParameter("leadPhoEtMin");
  string subleadPhoEtMinS  = cfg.getParameter("subleadPhoEtMin");
  doPtOverM = !(cfg.getParameter("doPtOverM").compare("no")==0);
  string PtOverMLeadS =    cfg.getParameter("PtOverMLead");
  string PtOverMSubLeadS  = cfg.getParameter("PtOverMSubLead");

  leadPhoEtMin = 33.;
  subleadPhoEtMin = 25.;
  PtOverMLead = 40./120.;
  PtOverMSubLead = 30./120.;
  if(leadPhoEtMinS.compare("") != 0) leadPhoEtMin = atof(leadPhoEtMinS.c_str());
  if(subleadPhoEtMinS.compare("") != 0) subleadPhoEtMin = atof(subleadPhoEtMinS.c_str());
  if(PtOverMLeadS.compare("") != 0) PtOverMLead = atof(PtOverMLeadS.c_str());
  if(PtOverMSubLeadS.compare("") != 0) PtOverMSubLead = atof(PtOverMSubLeadS.c_str());

  cout << "Kinematic Selections: " << endl 
       << "pt Lead:   " << leadPhoEtMin << "   pt SubLead:   " << subleadPhoEtMin << endl
       << "Do Pt Over M: " << doPtOverM << endl
       << "pt/m Lead: " << PtOverMLead <<  "   pt/m SubLead: " << PtOverMSubLead << endl;

  PhotonID = new HggPhotonID();
  PhotonID->setConfig(configFile);
  PhotonID->Init();

}

//
// Generic Initialization
//

int HggSelector::init(){
  clear();

  //setup MVA inputs histograms
  MVAInputs["massRes"] =  new TH1F("massRes","",200,0.,0.2);
  MVAInputs["massResWrongVtx"] = new TH1F("massResWrongVtx","",200,0.,0.2);
  MVAInputs["vtxProb"] = new TH1F("vtxProb","",100,0,1);
  MVAInputs["p1_EtByM"] = new TH1F("p1_EtByM","",200,0,1);
  MVAInputs["p2_EtByM"] = new TH1F("p2_EtByM","",200,0,1);
  MVAInputs["p1_Eta"] = new TH1F("p1_Eta","",120,-3,3);
  MVAInputs["p2_Eta"] = new TH1F("p2_Eta","",120,-3,3);
  MVAInputs["CosDPhi"] = new TH1F("cosDPhi","",200,-1,1);
  MVAInputs["p1_idMVA"] = new TH1F("p1_idMVA","",200,-1,1);
  MVAInputs["p2_idMVA"] = new TH1F("p2_idMVA","",200,-1,1);

  MjjDists["BeforeVBF"] = new TH1F("MjjBeforeVBF","",200,0,2000);
  MjjDists["dEta"] = new TH1F("dEta","",120,0,12);
  MjjDists["AfterDEta"] = new TH1F("MjjAfterDEta","",200,0,2000);
  MjjDists["Z"] = new TH1F("Z","",200,-20,20);
  MjjDists["AfterZ"] = new TH1F("MjjAfterZ","",200,0,2000);
  MjjDists["dPhi_jj_gg"] = new TH1F("dPhi_jj_gg","",32,0,3.2);
  MjjDists["Final"] = new TH1F("Final","",200,0,2000);

  this->setupTMVA();
  return 0;
}

void HggSelector::processEntry(Long64_t iEntry){
  processOnce();
  processEntry(kMVA);
  processEntry(kPFCiC);
  processEntry(kCiC);
}

void HggSelector::processOnce(){
  this->clear();
  if(!isData_) this->fillGenInfo();

  if(debugSelector) cout << "Setting Photon ID Vertices " << endl;
  PhotonID->setVertices(nVtx,vtxX,vtxY,vtxZ);

  if(debugSelector) cout << "requiring triggers ... " << endl;
  trigger_ = this->requireTrigger();

  // give the photon ID maker the vertex collection for this event

  if(debugSelector) cout << "Starting Photon Loop:  " << nPho_ << endl;
  for(int iPho=0; iPho<nPho_;iPho++){ //redo the scaling and smearing if needed
    VecbosPho *pho = &(Photons_->at(iPho));
    if(doScale && isData_){
      std::pair<float,float> dE = scale->getDEoE(*pho,runNumber);
      pho->dEoE    = dE.first;
      pho->dEoEErr = 0;
      pho->scaledEnergy = pho->correctedEnergy*(pho->dEoE);
      pho->scaledEnergyError = pho->correctedEnergyError*((pho->dEoE+pho->dEoEErr));
    }
    if(!isData_){
      pho->scaledEnergy = pho->correctedEnergy;
      pho->scaledEnergyError = pho->correctedEnergyError;
      if(doSmear){
	std::pair<float,float> dE = smear->getDEoE(*pho,applyScaleSmear);
	pho->dEoE    = dE.first;
	pho->dEoEErr = dE.second;
      }
      smearPhoton(pho,0);	
    }
  }


  MET = pfMet;
  METPhi = pfMetPhi;
  lumiBlockOut = lumiBlock;
  runNumberOut = runNumber;
  evtNumberOut = evtNumber;
  evtWeight = pileupWeight;

  eeBadScFilterFlagOut         = eeBadScFilterFlag;
  hcalLaserEventFilterFlagOut  = hcalLaserEventFilterFlag;
  HBHENoiseFilterResultFlagOut = HBHENoiseFilterResultFlag;
  isNotDeadEcalClusterOut      = isNotDeadEcalCluster;
  trackerFailureFilterFlagOut  = trackerFailureFilterFlag;
  CSCHaloFilterFlagOut         = CSCHaloFilterFlag;
  drDeadOut                    = drDead; 
  drBoundaryOut                = drBoundary;
  ECALTPFilterFlagOut          = ECALTPFilterFlag;

  nVtxOut = nVtx;
  rhoOut  = rho;

}

void HggSelector::processEntry(HggSelector::SelectionType SelType){
  int index1=-1,index2=-1;
  OutputVars *varPtr=&(vars[SelType]);
    
  float mvaOut[3];
  std::pair<int,int> indices;
  
  if(nSigma>0 && !isData_){ //do the energy scale and energy resolution systematics
    for(int iSmear=-nSigma;iSmear<=nSigma;iSmear++){

      switch(SelType){
      case kMVA:
	indices = getBestPair(mvaOut,iSmear,0);
	diPhoMVASmear.push_back(mvaOut[0]);
	pho1MVASmear.push_back(mvaOut[1]);
	pho2MVASmear.push_back(mvaOut[2]);
	break;
      case kPFCiC:
	indices = getBestPairCiC(iSmear,0,true);
	break;
      case kCiC:
	indices = getBestPairCiC(iSmear,0,false);
	break;
      default:
	std::cout << "Apparently this selection type isn't implemented, thats not good!" << std::endl;
	assert(false);
      }
	
      if(indices.first == -1 || indices.second==-1){
	varPtr->mPairSmear_.push_back(-1);
      }else{
	varPtr->mPairSmear_.push_back(getMPair(indices.first,indices.second));
      }
    }
  }

  for(int iScale=-nSigma;iScale<=nSigma;iScale++){ //do the scaling systematic
    
    
    switch(SelType){
    case kMVA:
      indices = getBestPair(mvaOut,-999,iScale);
      diPhoMVASmear.push_back(mvaOut[0]);
      pho1MVASmear.push_back(mvaOut[1]);
      pho2MVASmear.push_back(mvaOut[2]);
      break;
    case kPFCiC:
      indices = getBestPairCiC(0,iScale,true);
      break;
    case kCiC:
      indices = getBestPairCiC(0,iScale,false);
      break;
    default:
      std::cout << "Apparently this selection type isn't implemented, thats not good!" << std::endl;
      assert(false);
    }
    
    if(indices.first == -1 || indices.second==-1){
      varPtr->mPairScale_.push_back(-1);
    }else{
      varPtr->mPairScale_.push_back(getMPair(indices.first,indices.second));
    }
  }
  
  //check with the default smearing
  switch(SelType){
  case kMVA:
    indices = getBestPair(mvaOut,0,0); // no scaling and default smearing
    varPtr->diPhoMVA_=mvaOut[0];
    break;
  case kPFCiC:
    indices = getBestPairCiC(0,0,true); // no scaling and default smearing
    break;
  case kCiC:
    indices = getBestPairCiC(0,0,false); // no scaling and default smearing
    break;
  default:
    std::cout << "Apparently this selection type isn't implemented, thats not good!" << std::endl;
    assert(false);      
  }  
  
  index1 = indices.first;
  index2 = indices.second;
  
  if(debugSelector) cout << "LOOP DONE" << endl;	  
  
  
  if(debugSelector) cout << "indices: " << index1 << "  " << index2 << endl;
  
  
  if(index1 > -1 && index2 > -1){
    //fill MVA variables
    int selectedVertex = getVertexIndex(index1,index2);
    if(debugSelector) cout << "Final Selection MVA: " << selectedVertex << endl;
    
    TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
    pho1_ = Photons_->at(index1);
    pho2_ = Photons_->at(index2);
    varPtr->OutPhotons_.push_back(getReducedData(&pho1_,vtxPos,selectedVertex));
    varPtr->OutPhotons_.push_back(getReducedData(&pho2_,vtxPos,selectedVertex));
    
    if(debugSelector) cout << "Photon Indices: " << pho1_.index << "  " << pho2_.index << endl;
    
    varPtr->mPair_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
    varPtr->mPairNoCorr_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
    varPtr->mPairRes_ = massRes->getMassResolutionEonly(&pho1_,&pho2_,vtxPos);
    varPtr->mPairResWrongVtx_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
    varPtr->diPhoVtx_ = selectedVertex;
    varPtr->diPhoVtxX_ = vtxX[selectedVertex];
    varPtr->diPhoVtxY_ = vtxY[selectedVertex];
    varPtr->diPhoVtxZ_ = vtxZ[selectedVertex];
    varPtr->vtxProb_ = this->getVertexProb(index1,index2);
    
    float jpt[2];
    varPtr->Mjj_  = this->getVBFMjj(&pho1_,&pho2_,vtxPos,jpt);
    varPtr->ptJet1_ = jpt[0];
    varPtr->ptJet2_ = jpt[1];
    
    TLorentzVector p1 = pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy);
    TLorentzVector p2 = pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy);
    if(p1.Pt() < p2.Pt()){
      TLorentzVector tmp = p1;
      p1=p2; p2=tmp;
    }
    TLorentzVector gg = p1+p2;
    TVector3 boost = -1*gg.BoostVector();
    p1.Boost(boost);
    varPtr->cosThetaLead_ = p1.Vect().Dot(gg.Vect())/p1.Vect().Mag()/gg.Vect().Mag();
    
    if(SelType==kMVA) varPtr->cat_ = getCategory(*varPtr);
    else if(SelType == kPFCiC || kCiC) varPtr->cat_ = getCategoryPFCiC(*varPtr);
    else assert(false);
    
    float mva1 = PhotonID->getIdMVA(&pho1_,nVtx,rho,selectedVertex);
    float mva2 = PhotonID->getIdMVA(&pho2_,nVtx,rho,selectedVertex);
    
    for(int i1=0;i1<3;i1++){ for(int i2=0;i2<3;i2++){
	float offset1 = 0.01*(i1-1);
	float offset2 = 0.01*(i2-1);
	varPtr->diPhoMVAShift_[3*i1+i2] =  getDiPhoMVA(index1,index2,mva1+offset1,mva2+offset2,false);
      } }    
  }else{
    varPtr->mPair_=-1;      
    varPtr->cat_=-1;
  }

  varPtr->nOutPhotons_ = varPtr->OutPhotons_.size();
}

bool HggSelector::preSelectPhotons(VecbosPho* pho1,VecbosPho* pho2,TVector3 vtxPos){
  //apply kinematic photon selection
  if(fabs(pho1->SC.eta) > 2.5 || fabs(pho2->SC.eta) > 2.5) return false;  //outside of tracker acceptance
  if(fabs(pho1->SC.eta) > 1.4442 && fabs(pho1->SC.eta) < 1.566) return false;
  if(fabs(pho2->SC.eta) > 1.4442 && fabs(pho2->SC.eta) < 1.566) return false; // veto gap photons

  if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt() < subleadPhoEtMin || pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt() < subleadPhoEtMin ) return false;
  if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt() < leadPhoEtMin && pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt() < leadPhoEtMin ) return false;

  if(doPtOverM){
    float M = (pho1->p4FromVtx(vtxPos,pho1->finalEnergy) + pho2->p4FromVtx(vtxPos,pho2->finalEnergy)).M();
    if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()/M < PtOverMSubLead || pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt()/M < PtOverMSubLead ) return false;
    if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()/M < PtOverMLead && pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt()/M < PtOverMLead ) return false;
    
  }

  return true;
}

float HggSelector::getMPair(int i1, int i2){
  int selectedVertex = getVertexIndex(i1,i2);
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
  
  VecbosPho pho1 = Photons_->at(i1);
  VecbosPho pho2 = Photons_->at(i2);
  
  return (pho1.p4FromVtx(vtxPos,pho1.finalEnergy) + pho2.p4FromVtx(vtxPos,pho2.finalEnergy)).M();
}

std::pair<int,int> HggSelector::getBestPairCiC(int smearShift,int scaleShift,bool usePF=true){
  std::pair<int,int> indices(-1,-1);
  //int bestCat=-1;
  //float bestMass=-99;
  float highestPtSum=0; // is this the best way to select the pairs?
    for(int iPho1=0; iPho1<nPho_;iPho1++){
      if(photonMatchedElectron[iPho1] && doElectronVeto) continue;
      for(int iPho2=iPho1; iPho2<nPho_;iPho2++){
	if(iPho1==iPho2) continue;
	if(photonMatchedElectron[iPho2] && doElectronVeto) continue;
	if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
	//scale/smear the energy of the photon
	VecbosPho* pho1 = &(Photons_->at(iPho1));
	VecbosPho* pho2 = &(Photons_->at(iPho2));
	int selVtxI = this->getVertexIndex(iPho1,iPho2);
	TVector3 vtxPos(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]);
	if(!this->preSelectPhotons(pho1,pho2,vtxPos)) continue;
	if(!isData_){
	  
	  //apply scale shift	  
	  pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
	  pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;

	  smearPhoton(pho1,smearShift);
	  smearPhoton(pho2,smearShift);
	}
	bool CiC1,CiC2;
	if(usePF){
	  CiC1 = PhotonID->getIdCiCPF(pho1,nVtx,rho,selVtxI);
	  CiC2 = PhotonID->getIdCiCPF(pho2,nVtx,rho,selVtxI);
	}else{
	  CiC1 = PhotonID->getIdCiC(pho1,nVtx,rho,selVtxI);
	  CiC2 = PhotonID->getIdCiC(pho2,nVtx,rho,selVtxI);
	}
	if(!CiC1 || !CiC2) continue;
	float thisPtSum = pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()
	  + pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt();	
	if(thisPtSum > highestPtSum){
	  highestPtSum = thisPtSum;
	  indices.first = iPho1;
	  indices.second = iPho2;
	}
      }// for(iPho2...
    }// for(iPho1...
    return indices;
}

std::pair<int,int> HggSelector::getBestPair(float* mvaOut, int smearShift,int scaleShift){
  std::pair<int,int> indices(-1,-1);
  float diPhoMVAMax=-99;
  float maxSumPt=-1;
  *mvaOut = -999.;
  for(int iPho1=0; iPho1<nPho_;iPho1++){
    if(photonMatchedElectron[iPho1] && doElectronVeto) continue;
    for(int iPho2=iPho1; iPho2<nPho_;iPho2++){
      if(iPho1==iPho2) continue;
      if(photonMatchedElectron[iPho2] && doElectronVeto) continue;
      if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
      //NON-PF block
      //scale/smear the energy of the photon
      VecbosPho* pho1 = &(Photons_->at(iPho1));
      VecbosPho* pho2 = &(Photons_->at(iPho2));
      int selVtxI = this->getVertexIndex(iPho1,iPho2);
      TVector3 vtxPos(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]);
      if(!this->preSelectPhotons(pho1,pho2,vtxPos)) continue;
      if(!isData_){
	
	//apply scale shift	  
	pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
	pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;
	smearPhoton(pho1,smearShift);
	smearPhoton(pho2,smearShift);
      }
      if(debugSelector) cout << "Getting Photon ID:" << endl;
      float mva1 = PhotonID->getIdMVA(pho1,nVtx,rho,selVtxI);
      if(debugSelector) cout << "Getting Photon ID:" << endl;
      float mva2 = PhotonID->getIdMVA(pho2,nVtx,rho,selVtxI);
      if(debugSelector) cout << "Getting Double Photon MVA" << endl;
      float diPhoMVA =  getDiPhoMVA(iPho1,iPho2,mva1,mva2,false);
      if(debugSelector) cout << "\t\t" << mva1 << "  " << mva2 << "  " << diPhoMVA << endl;
      float thisPtSum = pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()
	+ pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt();	

      if(diPhoMVA>=-1 &&  thisPtSum > maxSumPt){
	maxSumPt = thisPtSum;
	indices.first = iPho1;
	indices.second = iPho2;
	diPhoMVAMax = diPhoMVA;
	*mvaOut = diPhoMVA;
      }
    }// for(iPho2...
  }// for(iPho1...
  
    return indices;
}

void HggSelector::smearPhoton(VecbosPho* pho,int smearShift){
  TRandom3 rng(0);
  float smear = pho->dEoE + smearShift*pho->dEoEErr;
  if(smear < 0) smear = 0;
  if(smearShift<-100) smear = 0;
  float rand=0;
  if(smear > 0) rand = rng.Gaus(0,smear);
  if(rand < -1) rand=-1;
  if(rand > 1E3) rand = 1E3;
  pho->finalEnergy = pho->scaledEnergy*(1+rand);
  pho->finalEnergyError = pho->scaledEnergyError;//*(1+rand);
}

void HggSelector::setDefaults(){
  clear();
}
void HggSelector::clear(){
  vars[kMVA].clear();
  vars[kPFCiC].clear();
  vars[kCiC].clear();

  pho1MVAScale.clear();
  pho2MVAScale.clear();
  diPhoMVAScale.clear();
  pho1MVASmear.clear();
  pho2MVASmear.clear();
  diPhoMVASmear.clear();
}

void HggSelector::fillGenInfo(){
  if(debugSelector) std::cout << "Setting Gen Higgs Info" << std::endl;
  if(nGenHiggs > 0){
    genHiggsPt = GenHiggs->front().pt;
    genHiggsVx = GenHiggs->front().Vx;
    genHiggsVy = GenHiggs->front().Vy;
    genHiggsVz = GenHiggs->front().Vz;
  }else{
    genHiggsPt = 0;
    genHiggsVx = 0;
    genHiggsVy = 0;
    genHiggsVz = 0;
  }

  if(debugSelector) std::cout << "Clearing Output Variables" << std::endl;
  int selGenPho=0;
  GenCollection::const_iterator genPho;
  if(GenPhotons == 0){
    std::cout << "WARNING: GenPhotons Collection is EMPTY!" << std::endl;
    return;
  }
  
  for(int i=0;i<nGenPho;i++){
    VecbosGen* genPho = &(GenPhotons->at(i));
    if(genPho->idMother != 25) continue; // not a higgs
    if(selGenPho==0){
      etaGenPho1 = genPho->eta;
      phiGenPho1 = genPho->phi;
      ptGenPho1 =  genPho->pt;
      energyGenPho1 = genPho->energy;
    }
    else if(selGenPho==1){
      etaGenPho2 = genPho->eta;
      phiGenPho2 = genPho->phi;
      ptGenPho2 =  genPho->pt;
      energyGenPho2 = genPho->energy;
    }else{
      cout << "More than 2 photons from a higgs decay!!" << endl;
    }
    selGenPho++;
  }
  
}

ReducedPhotonData HggSelector::getReducedData(VecbosPho* pho,TVector3 selVtx,int selVtxI){
  if(debugSelector) cout << "Filling Reduced Data: " << flush; 
  if(debugSelector) cout << pho->index << "  " << selVtx.Z() << "  " <<selVtxI <<endl;

  ReducedPhotonData data;
  TLorentzVector p4 = pho->p4FromVtx(selVtx,pho->finalEnergy,false);
  TLorentzVector p4NoCorr = pho->p4FromVtx(selVtx,pho->energy,false);
  TLorentzVector p4Gen;
  if(pho->genMatch.index!=-1) {
    p4Gen.SetPtEtaPhiE(pho->genMatch.pt,pho->genMatch.eta,pho->genMatch.phi,pho->genMatch.energy);
    data.pt_Gen = p4Gen.Pt(); data.eta_Gen = p4Gen.Eta(); data.phi_Gen = p4Gen.Phi(); data.E_Gen = p4Gen.E();
  }
  else{
    p4Gen.SetPtEtaPhiE(0.,0.,0.,0.);
    data.pt_Gen = 0.; data.eta_Gen = 0.; data.phi_Gen = 0.; data.E_Gen = 0.;  
  }
  data.pt = p4.Pt(); data.eta = p4.Eta(); data.phi = p4.Phi(); data.E =p4.E(); data.EError = pho->finalEnergyError;
  data.EErrorSmeared = massRes->getResolution(pho);
  data.index = pho->index;
  data.etaSC = pho->SC.eta;
  data.r9 = pho->SC.r9;
  data.passPFCiC = PhotonID->getIdCiCPF(pho,nVtx,rho,selVtxI); 
  data.category = (data.r9 < 0.94)+2*(fabs(data.etaSC) > 1.48); 
  data.idMVA = PhotonID->getIdMVA(pho,nVtx,rho,selVtxI);  
  data.mother = pho->genMatch.idMother; 
  PhotonID->fillIsoVariables(pho,&data,nVtx,rho,selVtxI);
  if(debugSelector) cout << "DONE Filling Reduced Data" <<endl;
  return data;
}

void HggSelector::setupTMVA(){
  diPhotonMVA = new TMVA::Reader( "!Color:!Silent" );; 

  // diphoton                
  diPhotonMVA->AddVariable("masserrsmeared/mass",&smearedMassErrByMass); 
  diPhotonMVA->AddVariable("masserrsmearedwrongvtx/mass",&smearedMassErrByMassWrongVtx);         
  diPhotonMVA->AddVariable("vtxprob",&vtxprob);             
  diPhotonMVA->AddVariable("ph1.pt/mass",&pho1PtByMass);         
  diPhotonMVA->AddVariable("ph2.pt/mass",&pho2PtByMass);         
  diPhotonMVA->AddVariable("ph1.eta",&pho1Eta);             
  diPhotonMVA->AddVariable("ph2.eta",&pho2Eta);             
  diPhotonMVA->AddVariable("TMath::Cos(ph1.phi-ph2.phi)",&cosDPhi);            
  diPhotonMVA->AddVariable("ph1.idmva",&pho1IdMVA);
  diPhotonMVA->AddVariable("ph2.idmva",&pho2IdMVA);

  //book MVAs:
  diPhotonMVA->BookMVA(  methodName_diPho, weightFile_diPho);
}

#define debugMjj 0
float HggSelector::getVBFMjj(VecbosPho* pho1, VecbosPho* pho2,TVector3 SelVtx,float *jetPts){
  TLorentzVector p1 = pho1->p4FromVtx(SelVtx,pho1->finalEnergy);
  TLorentzVector p2 = pho2->p4FromVtx(SelVtx,pho2->finalEnergy);
  jetPts[0]=0.; jetPts[1]=0.;
  if( max(p1.Pt(),p2.Pt())/(p1+p2).M() < 60./120. ) return 0;
  if( min(p1.Pt(),p2.Pt()) < 25. ) return 0;
  //if( (p1+p2).M() < 120. ) return 0; 
  if(debugMjj) cout << "NJets: " << nJet_ << endl;
  if(nJet_<2) return 0;
  int i1=-1,i2=-1;
  float maxSumPt=0;
  VecbosJet *jet1=0, *jet2=0;
  std::vector<VecbosJet>::iterator j1It;
  std::vector<VecbosJet>::iterator j2It;
  for(j1It = Jets_->begin(); j1It != Jets_->end(); j1It++){
    if(j1It->pt < 20.) continue;
    if(!passJetID(&*j1It)) continue;
    if(DeltaR(j1It->eta,p1.Eta(),j1It->phi,p1.Phi()) < 0.5 ) continue;
    if(DeltaR(j1It->eta,p2.Eta(),j1It->phi,p2.Phi()) < 0.5 ) continue;
    for(j2It = j1It+1; j2It != Jets_->end(); j2It++){
      if(j2It->pt < 20.) continue;
      if(!passJetID(&*j2It)) continue;
      if(DeltaR(j2It->eta,p1.Eta(),j2It->phi,p1.Phi()) < 0.5 ) continue;
      if(DeltaR(j2It->eta,p2.Eta(),j2It->phi,p2.Phi()) < 0.5 ) continue;

      if(j1It->pt < 30. && j2It->pt < 30.) continue;
      if(j1It->pt + j2It->pt > maxSumPt){
	maxSumPt = j1It->pt + j2It->pt;
	jet1 = &*j1It; jet2 = &*j2It;
      }
    }
  }

  if(jet1==0 || jet2==0) return 0;
  TLorentzVector j1 = jet1->getP4();
  TLorentzVector j2 = jet2->getP4();
  MjjDists["BeforeVBF"]->Fill((j1+j2).M());
  float dEtaJ = fabs(j1.Eta()-j2.Eta());
  if(debugMjj) cout << "dEtaJ: " << dEtaJ << endl;
  MjjDists["dEta"]->Fill(dEtaJ);
  if(dEtaJ < 3.) return 0;
  MjjDists["AfterDEta"]->Fill( (j1+j2).M() );
  TLorentzVector ggSystem = p1+p2;
  float Z = ggSystem.Eta() - (j1.Eta()+j2.Eta())/2.;
  if(debugMjj) cout << "Z: " << Z << endl;
  MjjDists["Z"]->Fill(Z);
  if(fabs(Z)>2.5) return 0;
  MjjDists["AfterZ"]->Fill( (j1+j2).M() );
  TLorentzVector jjSystem = j1+j2;
  if(debugMjj) cout << "dPhi jj gg: " << fabs(jjSystem.DeltaPhi(ggSystem)) << endl;
  MjjDists["dPhi_jj_gg"]->Fill( fabs(jjSystem.DeltaPhi(ggSystem)) );
  if( fabs(jjSystem.DeltaPhi(ggSystem)) < 2.6 ) return 0;
  MjjDists["Final"]->Fill( (j1+j2).M() );
  jetPts[0] = j1.Pt();
  jetPts[1] = j2.Pt();
  if(debugMjj) cout << "jj Mass: " << jjSystem.M() <<endl;
  return jjSystem.M();
}

bool HggSelector::passJetID(VecbosJet* jet){
  const int nJetCat=4;
  const float maxJetEta[nJetCat] = {2.5,2.75,3,4.7};
  const float betaStarSlope[nJetCat] = {0.2,0.3,999,999};
  const float rmsCut[nJetCat] = {0.06,0.05,0.05,0.055};

  int JetCat=-1;
  for(int i=0;i<nJetCat;i++){
    if(jet->eta < maxJetEta[i]){
      JetCat = i; break;
    }
  }
  if(JetCat == -1) return false;

  if(jet->betaStarClassicIdMVA > betaStarSlope[JetCat]*TMath::Log(nVtx-0.64)) return false; //warning, I think this is WRONG! but its how MIT does it...
  if(jet->rmsCandsHand > rmsCut[JetCat]) return false;
  return true;
}

int HggSelector::getCategory(OutputVars & vars){
  if(vars.diPhoMVA_ < -0.05) return -1;
  if(vars.Mjj_ >=500 && vars.ptJet1_ >= 30. && vars.ptJet2_ >= 30.) return 4;
  if(vars.Mjj_ >=250 && (vars.ptJet1_ >= 30. || vars.ptJet2_ >= 30.)) return 5;
  
  if(vars.diPhoMVA_ > 0.88) return 0;
  if(vars.diPhoMVA_ > 0.71) return 1;
  if(vars.diPhoMVA_ > 0.50) return 2;
  return 3;
}

int HggSelector::getCategoryPFCiC(OutputVars & vars){
  if(vars.Mjj_ >=500 && vars.ptJet1_ >= 30. && vars.ptJet2_ >= 30.) return 4;
  if(vars.Mjj_ >=250 && (vars.ptJet1_ >= 30. || vars.ptJet2_ >= 30.)) return 5;

  bool EBEB = (fabs(pho1_.SC.eta) < 1.48) && (fabs(pho2_.SC.eta) < 1.48);
  bool R9R9 = (pho1_.SC.r9 > 0.94) && (pho2_.SC.r9 > 0.94);
  
  if(R9R9 && EBEB) return 0;
  if(!R9R9 && EBEB) return 1;
  if(R9R9 && !EBEB) return 2;
  if(!R9R9 && !EBEB) return 3;
  return -1;
}


int HggSelector::getVertexIndex(int indexPho1,int indexPho2){
  if(forceVtxZero) return 0;
  //get the vertex selected by the di photon vertex MVA
  int selectedVertex = -1;
  if(indexPho1 > nPho_ || indexPho2 > nPho_) return -1;
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  if(debugSelector) cout << "Original Indices:" << origIndex1 << "  " << origIndex2 << endl;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      selectedVertex = ggVerticesVertexIndex->at(i).first;
      break;
    }
  }
  if(debugSelector) cout << "Returning Vertex Index: " << selectedVertex <<endl;
  return selectedVertex;
}
float HggSelector::getVertexMVA(int indexPho1,int indexPho2){
  //get the vertex selected by the di photon vertex MVA
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  float selectedVertex = -1;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      selectedVertex = ggVerticesVertexIndex->at(i).second;
      break;
    }
  }
  return selectedVertex;
}

float HggSelector::getVertexProb(int indexPho1,int indexPho2){
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  float evtmva = -1e6;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      evtmva = ggVerticesPerEvtMVA->at(i);
      break;
    }
  }
  return vtxprob=1.-0.49*(evtmva+1.0);
  
}

float HggSelector::getDiPhoMVA(int indexPho1, int indexPho2, float mva1, float mva2, bool usePFSC){
  if(debugSelector) cout << "getDiPhoMVA" <<endl;  
  int selectedVertex = getVertexIndex(indexPho1,indexPho2);
  if(selectedVertex<0 || selectedVertex >= nVtx){
    cout << "WARNING: Photons " << indexPho1 << " and " << indexPho2 
	 << " have no selected vertex!" << endl
	 << "Skipping this pair" << endl;
    return -9999.;
  }
  if(debugSelector) cout << mva1 << "  " << mva2 << endl;
  if(mva1 < -999 || mva2 < -999) return -9999.;
  if(debugSelector) cout << "getting photons" <<endl;  
  VecbosPho pho1 = Photons_->at(indexPho1);
  VecbosPho pho2 = Photons_->at(indexPho2);
  if(debugSelector) cout << "selected Vertex: " << selectedVertex <<endl;  

  if(debugSelector){
    std::cout << "pho1: " << pho1.finalEnergy << "  " << pho1.scaledEnergy << "  " << pho1.dEoE << "  " << pho1.eta << "  " << std::endl;
    std::cout << "pho2: " << pho2.finalEnergy << "  " << pho2.scaledEnergy << "  " << pho2.dEoE << "  " << pho2.eta << "  " << std::endl;
  }
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);

  TLorentzVector p1 = pho1.p4FromVtx(vtxPos,pho1.finalEnergy);
  TLorentzVector p2 = pho2.p4FromVtx(vtxPos,pho2.finalEnergy);

  VecbosPho* phoLead = (p1.Et() > p2.Et() ? &pho1 : & pho2); //select the leading photon
  VecbosPho* phoSubLead = (p1.Et() > p2.Et() ? &pho2 : & pho1); //select the leading photon

  VecbosSC *scLead;
  VecbosSC *scSubLead;
  if(usePFSC){
    scLead = &(phoLead->PFSC);
    scSubLead = &(phoSubLead->PFSC);
  }else{
    scLead = &(phoLead->SC);
    scSubLead = &(phoSubLead->SC);
  }

  float mvaLead = (p1.Et() > p2.Et() ? mva1 : mva2);
  float mvaSubLead = (p1.Et() > p2.Et() ? mva2 : mva1);

  float mPair = (p1+p2).M();
  if(debugSelector) std::cout << "mPair: " << mPair << std::endl;
  //fill variables
  smearedMassErrByMass = massRes->getMassResolutionEonly(&pho1,&pho2,vtxPos)/mPair;
  smearedMassErrByMassWrongVtx = massRes->getMassResolution(&pho1,&pho2,vtxPos,true)/mPair;
  if(debugSelector) cout << "Mass Erro resolved. Selected Vertex: " << selectedVertex << endl;

  vtxprob= getVertexProb(indexPho1,indexPho2);


  if(debugSelector) cout << "vtx prob: " << vtxprob << endl;
  pho1PtByMass = max(p1.Et(),p2.Et())/mPair;
  pho2PtByMass = min(p1.Et(),p2.Et())/mPair;
  if(debugSelector) cout << "pho1/m: " << pho1PtByMass << endl;
  if(debugSelector) cout << "pho2/m: " << pho2PtByMass << endl;
  pho1Eta = (p1.Et() > p2.Et() ? p1.Eta(): p2.Eta());
  pho2Eta = (p1.Et() > p2.Et() ? p2.Eta(): p1.Eta());
  if(debugSelector) cout << "pho1 eta: " << pho1Eta << endl;
  if(debugSelector) cout << "pho2 eta: " << pho2Eta << endl;
  cosDPhi = TMath::Cos(DeltaPhi(scLead->phi,scSubLead->phi));
  if(debugSelector) cout << "cos dPhi: " << cosDPhi << endl;
  pho1IdMVA = mvaLead;
  pho2IdMVA = mvaSubLead;

  MVAInputs["massRes"]->Fill(smearedMassErrByMass);
  MVAInputs["massResWrongVtx"]->Fill(smearedMassErrByMassWrongVtx);
  MVAInputs["vtxProb"]->Fill(vtxprob);
  MVAInputs["p1_EtByM"]->Fill(pho1PtByMass);
  MVAInputs["p2_EtByM"]->Fill(pho2PtByMass);
  MVAInputs["p1_Eta"]->Fill(pho1Eta);
  MVAInputs["p2_Eta"]->Fill(pho2Eta);
  MVAInputs["CosDPhi"]->Fill(cosDPhi);
  MVAInputs["p1_idMVA"]->Fill(pho1IdMVA);
  MVAInputs["p2_idMVA"]->Fill(pho2IdMVA);

  return diPhotonMVA->EvaluateMVA(methodName_diPho);
}

void HggSelector::setupOutputTree(){
  outTree = new TTree("HggOutput","");
  outTree->Branch("trigger",&trigger_,"trigger/I");


  //MVA Selection
  outTree->Branch("mPair",&vars[0].mPair_,"mPair/F");
  outTree->Branch("mPairNoCorr",&vars[0].mPairNoCorr_,"mPairNoCorr/F");
  outTree->Branch("mPairRes",&vars[0].mPairRes_,"mPairRes/F");
  outTree->Branch("mPairResWrongVtx",&vars[0].mPairResWrongVtx_,"mPairResWrongVtx/F");
  outTree->Branch("diPhotonMVA",&vars[0].diPhoMVA_,"diPhotonMVA/F");
  outTree->Branch("diPhotonMVAShift",vars[0].diPhoMVAShift_,"diPhotonMVAShift[9]/F");
  outTree->Branch("diPhotonVtx",&vars[0].diPhoVtx_,"diPhotonVtx/I");
  outTree->Branch("diPhotonVtxX",&vars[0].diPhoVtxX_,"diPhotonVtxX/F");
  outTree->Branch("diPhotonVtxY",&vars[0].diPhoVtxY_,"diPhotonVtxY/F");
  outTree->Branch("diPhotonVtxZ",&vars[0].diPhoVtxZ_,"diPhotonVtxZ/F");
  outTree->Branch("vtxProb",&vars[0].vtxProb_,"vtxProb/F");
  outTree->Branch("Mjj",&vars[0].Mjj_,"Mjj");
  outTree->Branch("ptJet1",&vars[0].ptJet1_,"ptJet1");
  outTree->Branch("ptJet2",&vars[0].ptJet2_,"ptJet2");
  outTree->Branch("cosThetaLead",&vars[0].cosThetaLead_,"cosThetaLead/F");
  outTree->Branch("cat",&vars[0].cat_,"cat/I");
  
  outTree->Branch("mPairScale",&vars[0].mPairScale_);
  outTree->Branch("pho1MVAScale",&pho1MVAScale);
  outTree->Branch("pho2MVAScale",&pho2MVAScale);
  outTree->Branch("diPhoMVAScale",&diPhoMVAScale);
  outTree->Branch("mPairSmear",&vars[0].mPairSmear_);
  outTree->Branch("pho1MVASmear",&pho1MVASmear);
  outTree->Branch("pho2MVASmear",&pho2MVASmear);
  outTree->Branch("diPhoMVASmear",&diPhoMVASmear);

  outTree->Branch("nPhoton",&vars[0].nOutPhotons_,"nPhoton/I");
  outTree->Branch("Photon",&vars[0].OutPhotons_);


  outTree->Branch("mPairPFCiC",&vars[1].mPair_,"mPairPFCiC/F");
  outTree->Branch("mPairNoCorrPFCiC",&vars[1].mPairNoCorr_,"mPairNoCorrPFCiC/F");
  outTree->Branch("mPairResPFCiC",&vars[1].mPairRes_,"mPairResPFCiC/F");
  outTree->Branch("mPairResWrongVtxPFCiC",&vars[1].mPairResWrongVtx_,"mPairResWrongVtxPFCiC/F");
  outTree->Branch("diPhotonVtxPFCiC",&vars[1].diPhoVtx_,"diPhotonVtxPFCiC/I");
  outTree->Branch("diPhotonVtxXPFCiC",&vars[1].diPhoVtxX_,"diPhotonVtxXPFCiC/F");
  outTree->Branch("diPhotonVtxYPFCiC",&vars[1].diPhoVtxY_,"diPhotonVtxYPFCiC/F");
  outTree->Branch("diPhotonVtxZPFCiC",&vars[1].diPhoVtxZ_,"diPhotonVtxZPFCiC/F");
  outTree->Branch("vtxProbPFCiC",&vars[1].vtxProb_,"vtxProbPFCiC/F");
  outTree->Branch("MjjPFCiC",&vars[1].Mjj_,"MjjPFCiC");
  outTree->Branch("ptJet1PFCiC",&vars[1].ptJet1_,"ptJet1PFCiC");
  outTree->Branch("ptJet2PFCiC",&vars[1].ptJet2_,"ptJet2PFCiC");
  outTree->Branch("cosThetaLeadPFCiC",&vars[1].cosThetaLead_,"cosThetaLeadPFCiC/F");
  outTree->Branch("catPFCiC",&vars[1].cat_,"catPFCiC/I");

  outTree->Branch("mPairScalePFCiC",&vars[1].mPairScale_);
  outTree->Branch("mPairSmearPFCiC",&vars[1].mPairSmear_);

  outTree->Branch("nPhotonPFCiC",&vars[1].nOutPhotons_,"nPhotonPFCiC/I");
  outTree->Branch("PhotonPFCiC",&vars[1].OutPhotons_);

  //CiC Selection
  outTree->Branch("mPairCiC",&vars[2].mPair_,"mPairCiC/F");
  outTree->Branch("mPairNoCorrCiC",&vars[2].mPairNoCorr_,"mPairNoCorrCiC/F");
  outTree->Branch("mPairResCiC",&vars[2].mPairRes_,"mPairResCiC/F");
  outTree->Branch("mPairResWrongVtxCiC",&vars[2].mPairResWrongVtx_,"mPairResWrongVtxCiC/F");
  outTree->Branch("diPhotonVtxCiC",&vars[2].diPhoVtx_,"diPhotonVtxCiC/I");
  outTree->Branch("diPhotonVtxXCiC",&vars[2].diPhoVtxX_,"diPhotonVtxXCiC/F");
  outTree->Branch("diPhotonVtxYCiC",&vars[2].diPhoVtxY_,"diPhotonVtxYCiC/F");
  outTree->Branch("diPhotonVtxZCiC",&vars[2].diPhoVtxZ_,"diPhotonVtxZCiC/F");
  outTree->Branch("vtxProbCiC",&vars[2].vtxProb_,"vtxProbCiC/F");
  outTree->Branch("MjjCiC",&vars[2].Mjj_,"MjjCiC");
  outTree->Branch("ptJet1CiC",&vars[2].ptJet1_,"ptJet1CiC");
  outTree->Branch("ptJet2CiC",&vars[2].ptJet2_,"ptJet2CiC");
  outTree->Branch("cosThetaLeadCiC",&vars[2].cosThetaLead_,"cosThetaLeadCiC/F");


  outTree->Branch("mPairScaleCiC",&vars[2].mPairScale_);
  outTree->Branch("mPairSmearCiC",&vars[2].mPairSmear_);

  outTree->Branch("nPhotonCiC",&vars[2].nOutPhotons_,"nPhotonCiC/I");
  outTree->Branch("PhotonCiC",&vars[2].OutPhotons_);



  outTree->Branch("eeBadScFilterFlag",&eeBadScFilterFlagOut);
  outTree->Branch("hcalLaserEventFilterFlag",&hcalLaserEventFilterFlagOut);
  outTree->Branch("HBHENoiseFilterResultFlag",&HBHENoiseFilterResultFlagOut);
  outTree->Branch("isNotDeadEcalCluster",&isNotDeadEcalClusterOut);
  outTree->Branch("trackerFailureFilterFlag",&trackerFailureFilterFlagOut);
  outTree->Branch("CSCHaloFilterFlag",&CSCHaloFilterFlagOut);
  outTree->Branch("drDead",&drDeadOut);
  outTree->Branch("drBoundary",&drBoundaryOut);
  outTree->Branch("ECALTPFilterFlag",&ECALTPFilterFlagOut);


  outTree->Branch("MET",&MET);
  outTree->Branch("METPhi",&METPhi);

  outTree->Branch("genHiggsPt",&genHiggsPt);
  outTree->Branch("genHiggsVx",&genHiggsVx);
  outTree->Branch("genHiggsVy",&genHiggsVy);
  outTree->Branch("genHiggsVz",&genHiggsVz);
  
  outTree->Branch("evtWeight",&evtWeight);

  outTree->Branch("ptGenPho1",&ptGenPho1);
  outTree->Branch("etaGenPho1",&etaGenPho1);
  outTree->Branch("phiGenPho1",&phiGenPho1);
  outTree->Branch("energyGenPho1",&energyGenPho1);

  outTree->Branch("ptGenPho2",&ptGenPho2);
  outTree->Branch("etaGenPho2",&etaGenPho2);
  outTree->Branch("phiGenPho2",&phiGenPho2);
  outTree->Branch("energyGenPho2",&energyGenPho2);


  outTree->Branch("nPU",&nPU_);
  outTree->Branch("nVtx",&nVtxOut);
  outTree->Branch("rho",&rhoOut);

  outTree->Branch("runNumber",&runNumberOut);
  outTree->Branch("evtNumber",&evtNumberOut);
  outTree->Branch("lumiBlock",&lumiBlockOut);
  }
