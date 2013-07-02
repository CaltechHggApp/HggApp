#include "MakeSpinWorkspace.h"
#include <fstream>

#include "assert.h"

using namespace std;
#include "selectionMaps.C"

#define debug_workspaces 0

MakeSpinWorkspace::MakeSpinWorkspace(TString outputFileName):
  mixer(0),
  fileKFactor(0),
  fileRescaleFactor(0),
  isGlobe(0),
  lumi(-1)
{

  optimization=false;
  //setup the defaults
  takeCatFromTree=false;
  twoEBcats=false;
  vetoInnerEE=false;

  requireCiC=true;
  tightPt=false;
  selectionMap=0;

  useUncorrMass=false;

  runLow=0;
  runHigh= 9999999;

  useR9=false;
  nCat=4;

  setMassRange(110.,170.);
  setPtCuts(32.,24.);

  setUseHelicityFrame(false);
  setUseAbsCosTheta(false);

  EfficiencyCorrectionFile_Data="";
  EfficiencyCorrectionFile_MC="";
  filenameKFactor="";
  filenameRescaleFactor="";

  //open output file and create the new workspace
  outputFile = TFile::Open(outputFileName,"RECREATE");
  ws = new RooWorkspace();
  ws->SetName("cms_hgg_spin_workspace");  
  labels = new RooCategory("labels","labels");
  
  /*
  chargedIso = {3.8,2.5,3.1,2.2};
  goodIsoSum = {6.0,4.7,5.6,3.6};
  badIsoSum  = {10.,6.5,5.6,4.4};
  */

  chargedIso[0] = 3.8; chargedIso[1] = 2.5; chargedIso[2] = 3.1; chargedIso[3] = 2.2;
  goodIsoSum[0] = 6.0; goodIsoSum[1] = 4.7; goodIsoSum[2] = 5.6; goodIsoSum[3] = 3.6;
  badIsoSum [0] = 10.; badIsoSum [1] = 6.5; badIsoSum [2] = 5.6; badIsoSum [3] = 4.4;
}

MakeSpinWorkspace::~MakeSpinWorkspace(){
  if(mixer) delete mixer;
}


void MakeSpinWorkspace::getSelectionMap(int map,bool isData){
  // get selection map
  // needs to be cleaned up, but useful for testing atm

  selectionMaps.clear();
  switch(map){
    /* // OLD MAPS
  case 0:
    selectionMaps.push_back(getSelectionMap0());
  case 1:
    selectionMaps.push_back(getSelectionMap1());
  case 2:
    selectionMaps.push_back(getSelectionMap2());
  case 3:
    selectionMaps.push_back(getSelectionMap3());
  case 4:
    selectionMaps.push_back(getSelectionMap4(isData));
  case 5:
    selectionMaps.push_back(getSelectionMap5(isData));
  case 6:
    selectionMaps.push_back(getSelectionMap6(isData));
  case 7:
    selectionMaps.push_back(getSelectionMap7(isData));
    */
  case 0:
    selectionMaps = getSelectionMap8();    
    break;
  case 1:
    selectionMaps = getSelectionMap9();    
    break;
  case 2:
    selectionMaps = getSelectionMap10();    
    break;
  case 3:
    selectionMaps = getSelectionMap11();    
    break;
  case 4:
    selectionMaps = getSelectionMap12();    
    break;
  case 5:
    selectionMaps = getSelectionMap13();    
    break;
  case 6:
    selectionMaps = getSelectionMap14();    
    break;
  }
  nCat = 2*selectionMaps.size();
  if(twoEBcats) nCat += selectionMaps.size();
}

int MakeSpinWorkspace::passSelection(TH2F* map,float sigEoE1,float eta1, float pt1,float sigEoE2,float eta2, float pt2){
  float cut1 = map->GetBinContent(map->FindFixBin(fabs(eta1),pt1)); 
  float cut2 = map->GetBinContent(map->FindFixBin(fabs(eta2),pt2)); 
  if(sqrt(sigEoE1*sigEoE1/cut1/cut1+sigEoE2*sigEoE2/cut2/cut2) <=sqrt(2)) return 0; //elliptical cut!
  return 1;
}

int MakeSpinWorkspace::passSelection(float r9){
  if(r9 > 0.94) return 0;
  return 1;
}
bool MakeSpinWorkspace::getBaselineSelection(HggOutputReader2* h,int maxI,int minI,float mass){
  if(requireCiC){
  if(h->nPhotonPFCiC < 2
     || h->diPhotonMVA<-1
     || mass < mMin
     || mass > mMax
     || h->PhotonPFCiC_pt[minI] < pt2Min
     || h->PhotonPFCiC_pt[maxI] < pt1Min) return false;
  }else{
    if(h->nPhoton < 2
       || h->diPhotonMVA<-1
       || mass < mMin
       || mass > mMax
       || h->Photon_pt[minI] < pt2Min
       || h->Photon_pt[maxI] < pt1Min) return false;
  }
  return true;
}

void MakeSpinWorkspace::AddToWorkspace(TString inputFile,TString tag, bool isData, int N, bool isList){
  //adds the data/MC to the workspace with the given tag
  std::cout << "Processing " <<tag << " : " << inputFile <<std::endl;
  TFile *f=0;
  TChain *tree=0;
  TString treeName = (isGlobe ? Form("spin_trees/%s",tag.Data()) : "HggOutput");
  if(isList){
    std::cout << "isList" <<std::endl;
    tree = getChainFromList(inputFile,treeName);
  }else{
    f = TFile::Open(inputFile);
    tree = (TChain*)f->Get(treeName);
  }
  if(tree==0) return;

  //define a reader class for the input tree
  HggOutputReader2 *h=0;
  GlobeReader *g=0;
  if(isGlobe) g = new GlobeReader(tree);
  else h = new HggOutputReader2(tree);
    

  if(!isGlobe) setupBranches(*h);
  
  TFile *efficiencyCorrection=0;
  if(isData && EfficiencyCorrectionFile_Data!=""){ //Use MC derived correction weights for the photons
    //as a function of eta/R9/phi
    std::cout << "Using Efficiency Correction " << EfficiencyCorrectionFile_Data <<std::endl;
    efficiencyCorrection = new TFile(EfficiencyCorrectionFile_Data);
  }
  if(!isData && EfficiencyCorrectionFile_MC!=""){ //Use MC derived correction weights for the photons
    //as a function of eta/R9/phi
    std::cout << "Using Efficiency Correction " << EfficiencyCorrectionFile_MC <<std::endl;
    efficiencyCorrection = new TFile(EfficiencyCorrectionFile_MC);
  }
  effMaps.clear();
  if(efficiencyCorrection){
    TList *keyList =  efficiencyCorrection->GetListOfKeys();
    for(int i=0;i<keyList->GetEntries();i++){
      effMaps.push_back( (TH3F*)efficiencyCorrection->Get(keyList->At(i)->GetName()) );
    }

  }
  if(filenameKFactor!=""){
    if(tag=="Hgg125") fileKFactor = new TFile(filenameKFactor);
    else fileKFactor=0;
  }

  if(filenameRescaleFactor!=""){
    fileRescaleFactor = new TFile(filenameRescaleFactor);
  }

  //define the variables for the workspace
  RooRealVar* mass   = new RooRealVar("mass",  "Mass [GeV]", mMin, mMax);
  RooRealVar* cosT   = new RooRealVar("cosT",  "cos(theta)", (useAbsCosTheta ? 0 : -1),1);
  RooRealVar* sige1  = new RooRealVar("sigEoE1","#sigma_{E}/E Lead Photon",0,0.1);
  RooRealVar* sige2  = new RooRealVar("sigEoE2","#sigma_{E}/E SubLead Photon",0,0.1);
  RooRealVar* evtW   = new RooRealVar("evtWeight","Event Weight",1,0,20);

  RooRealVar* eta1 = new RooRealVar("eta1","#eta Lead Photon",0,-3,3);
  RooRealVar* eta2 = new RooRealVar("eta2","#eta SubLead Photon",0,-3,3);

  RooRealVar* etaSC1 = new RooRealVar("etaSC1","SC #eta Lead Photon",0,-3,3);
  RooRealVar* etaSC2 = new RooRealVar("etaSC2","SC #eta SubLead Photon",0,-3,3);

  RooRealVar* phi1 = new RooRealVar("phi1","#phi Lead Photon",0,0,6.3);
  RooRealVar* phi2 = new RooRealVar("phi2","#phi SubLead Photon",0,0,6.3);

  RooRealVar* pt1 = new RooRealVar("pt1","p_{T} Lead Photon",0,0,2e3);
  RooRealVar* pt2 = new RooRealVar("pt2","p_{T} SubLead Photon",0,0,2e3);

  RooRealVar* r91 = new RooRealVar("r91","R_{9} Lead Photon",0,0,1.3);
  RooRealVar* r92 = new RooRealVar("r92","R_{9} SubLead Photon",0,0,1.3);

  RooRealVar* idMVA1 = new RooRealVar("idMVA1","ID MVA Lead Photon",0,-1,1.);
  RooRealVar* idMVA2 = new RooRealVar("idMVA2","ID MVA SubLead Photon",0,-1,1.);

  RooRealVar* diPhotonMVA = new RooRealVar("diPhotonMVA","DiPhoton MVA",0,-1,1);

  RooCategory *cat = new RooCategory("evtcat","evtcat");

  //create the RooArgSet to initialize the datasets
  RooArgSet set;
  set.add(*mass);
  set.add(*cosT);
  set.add(*evtW);

  if(optimization){
    set.add(*sige1);
    set.add(*sige2);
    
    set.add(*eta1); set.add(*eta2);
    set.add(*etaSC1); set.add(*etaSC2);
    set.add(*phi1); set.add(*phi2);
    set.add(*pt1); set.add(*pt2);
    set.add(*r91); set.add(*r92);
  }
  /*
  set.add(*idMVA1); set.add(*idMVA2);
  set.add(*diPhotonMVA);
  */
  RooRealVar *totEB = new RooRealVar(Form("%s_EB_totalEvents",tag.Data()),"",0,0,1e9);
  RooRealVar *totEE = new RooRealVar(Form("%s_EE_totalEvents",tag.Data()),"",0,0,1e9);
  RooRealVar *totGen = new RooRealVar(tag+"_Ngen","",N);

  std::map<int,RooDataSet*> dataMap;

  std::cout << "nCat: " << nCat <<std::endl;
  for(int j=0;j<nCat;j++){
    dataMap[j] = new RooDataSet(Form("%s_cat%d",tag.Data(),j),"",set);    
  }

  Long64_t iEntry=-1;
  cout << "Making DataSet" << endl;
  double nEB=0,nEE=0; // store the total number of events before any cuts or selection

  //need to correct for overweighting the signal MC for the PU (<w> != 1)
  long int Nentries=0;

  //begin main loop over tree
  while(true){
    int val=1;
    if(isGlobe) val=g->GetEntry(++iEntry);
    else val=h->GetEntry(++iEntry);
    
    if(val==0) break;
    if(debug_workspaces && iEntry>=20000) break;

    if( !(iEntry%10000) ){
      cout << "Processing " << iEntry << endl;

    }
    //determine  the leading and trailing photons
    int maxI,minI;
    if(isGlobe){
      maxI=0; minI=1;
    }else{
      if(requireCiC){
	maxI = (h->PhotonPFCiC_pt[1] > h->PhotonPFCiC_pt[0] ? 1:0);
	minI = (h->PhotonPFCiC_pt[1] > h->PhotonPFCiC_pt[0] ? 0:1);
      }else{
	maxI = (h->Photon_pt[1] > h->Photon_pt[0] ? 1:0);
	minI = (h->Photon_pt[1] > h->Photon_pt[0] ? 0:1);
      }
    }

    if(!isGlobe){ // globe trees don't have run number ...
      int run = h->runNumber;
      if(isData && (run < runLow || run > runHigh) ) continue;
    }
    float pho1_etaSC, pho2_etaSC;
    float se1,se2;
    float r91_f,r92_f;
    float pt1_f,pt2_f;
    float eta1_f,eta2_f;
    float phi1_f,phi2_f;
    if(isGlobe){
      pho1_etaSC = g->lead_calo_eta;
      pho2_etaSC = g->sublead_calo_eta;
      se1 = g->lead_sigmaE_nosmear/g->lead_E;
      se2 = g->sublead_sigmaE_nosmear/g->sublead_E;
      r91_f = g->lead_r9;
      r92_f = g->sublead_r9;
      pt1_f = TMath::Sqrt(TMath::Power(g->lead_px,2)+TMath::Power(g->lead_py,2));
      pt2_f = TMath::Sqrt(TMath::Power(g->sublead_px,2)+TMath::Power(g->sublead_py,2));
      eta1_f = g->lead_calo_eta;
      eta2_f = g->sublead_calo_eta;
      phi1_f = g->lead_calo_phi;
      phi2_f = g->sublead_calo_phi;
    }else{
      if(requireCiC){
	pho1_etaSC = h->PhotonPFCiC_etaSC[maxI];
	pho2_etaSC = h->PhotonPFCiC_etaSC[minI];
	se1 = h->PhotonPFCiC_EError[maxI]/h->PhotonPFCiC_E[maxI];
	se2 = h->PhotonPFCiC_EError[minI]/h->PhotonPFCiC_E[minI];
	r91_f = h->PhotonPFCiC_r9[maxI];
	r92_f = h->PhotonPFCiC_r9[minI];
	pt1_f = h->PhotonPFCiC_pt[maxI];
	pt2_f = h->PhotonPFCiC_pt[minI];
	eta1_f = h->PhotonPFCiC_etaSC[maxI];
	eta2_f = h->PhotonPFCiC_etaSC[minI];
	phi1_f = h->PhotonPFCiC_phi[maxI];
	phi2_f = h->PhotonPFCiC_phi[minI];
      }else{
	pho1_etaSC = h->Photon_etaSC[maxI];
	pho2_etaSC = h->Photon_etaSC[minI];
	se1 = h->Photon_EError[maxI]/h->Photon_E[maxI];
	se2 = h->Photon_EError[minI]/h->Photon_E[minI];
	r91_f = h->Photon_r9[maxI];
	r92_f = h->Photon_r9[minI];
	pt1_f = h->Photon_pt[maxI];
	pt2_f = h->Photon_pt[minI];
	eta1_f = h->Photon_etaSC[maxI];
	eta2_f = h->Photon_etaSC[minI];
	phi1_f = h->Photon_phi[maxI];
	phi2_f = h->Photon_phi[minI];
      }
    }
    float weight=1;
    if(!isData) weight = (isGlobe ? g->evweight : h->evtWeight*getEfficiency(*h,125));
    if(weight==0 && !isData && !isGlobe){
      std::cout << "0 weight event!    " << h->evtWeight << "    " << getEfficiency(*h,125) <<std::endl;
    }
    if(fabs(pho1_etaSC) < 1.48 && fabs(pho2_etaSC) < 1.48) nEB+=weight;
    else nEE+=weight;

    Nentries++;
    
    float m;
    if(isGlobe) m = g->higgs_mass;
    else{
      if(requireCiC) m = (useUncorrMass ? h->mPairNoCorrPFCiC : h->mPairPFCiC);
      else m = (useUncorrMass ? h->mPairNoCorr : h->mPair);

    }
    if(debug_workspaces && m!=-1) std::cout << "m: " << m <<std::endl;

    // apply selections
    if(isGlobe){ 
      if(m <mMin || m > mMax) continue;
    }else{// globe ntuples already have selections, just need mass cuts
      if(!getBaselineSelection(h,maxI,minI,m)) continue;
      if(tightPt && (pt1_f/m < 1./3. || pt2_f/m < 1./4.) ) continue;
    }
    if(debug_workspaces) std::cout <<  "passed selection" <<std::endl;

    int p1=nCat,p2=nCat;
    assert(minI != maxI);
    //determine the photon categories
    int cat=-1;
    if(vetoInnerEE){
      if(fabs(pho1_etaSC) > 1.48 && fabs(pho1_etaSC) < 1.70) continue;
      if(fabs(pho2_etaSC) > 1.48 && fabs(pho2_etaSC) < 1.70) continue;
    }
    if(takeCatFromTree && isGlobe){
      cat = g->category;
    }else if(useR9){
      p1 = passSelection(r91_f);
      p2 = passSelection(r92_f);
      if(p1 && p2) cat=0;
      else cat=1;
      if(twoEBcats) if( !(fabs(pho1_etaSC) < 0.80 && fabs(pho2_etaSC) < 0.80) ) cat+=2;
      if( !(fabs(pho1_etaSC) < 1.48 && fabs(pho2_etaSC) < 1.48) ) cat+=2;
    }else{
      if(!isGlobe) if(!passCiCIso(*h,0) || !passCiCIso(*h,1)) continue;// veto if either photon fails the iso

      int nLayer=nCat/(twoEBcats?3:2);
      for(int iCut=0;iCut<nLayer;iCut++){
	if(passSelection(selectionMaps.at(iCut),se1,fabs(eta1_f),pt1_f,se2,fabs(eta2_f),pt2_f) == 0){
	  cat=iCut;
	  break;
	}
      }
      if(twoEBcats) if( !(fabs(pho1_etaSC) < 0.80 && fabs(pho2_etaSC) < 0.80) ) cat+=nLayer;
      if( !(fabs(pho1_etaSC) < 1.48 && fabs(pho2_etaSC) < 1.48) ) cat+=nLayer;      
    }
    if(cat < 0 || cat >= nCat) continue;
    //set all the variables
    mass->setVal(m);
    sige1->setVal(se1);
    sige2->setVal(se2);

    eta1->setVal(eta1_f);
    eta2->setVal(eta2_f);
    etaSC1->setVal(pho1_etaSC);
    etaSC2->setVal(pho2_etaSC);
    phi1->setVal(phi1_f);
    phi2->setVal(phi2_f);
    pt1->setVal(pt1_f);
    pt2->setVal(pt2_f);
    r91->setVal(r91_f);
    r92->setVal(r92_f);

    float cosTheta;    

    if(isGlobe){
      if(useHelicityFrame) cosTheta = g->costheta_hx;
      else                 cosTheta = g->costheta_cs;
    }else{ // not glob
      if(useHelicityFrame) cosTheta = h->cosThetaLead;
      else                 cosTheta = calculateCosThetaCS(h);
    }
    if(useAbsCosTheta) cosTheta = fabs(cosTheta);
    cosT->setVal(cosTheta);

    //efficiencyCorrection=0
    float pho1EffWeight = getEffWeight(eta1_f,pt1_f,phi1_f,r91_f);
    float pho2EffWeight = getEffWeight(eta2_f,pt2_f,phi2_f,r92_f);

    //set the event weight
    if(!isData) evtW->setVal(weight*pho1EffWeight*pho2EffWeight);
    else evtW->setVal(1*pho1EffWeight*pho2EffWeight);
    if( !(iEntry%1000) ) cout <<  "\t\t\t" << evtW->getVal() << endl;


    //if(evtW->getVal() > 7.4) cout << iEntry << ": " << mass->getVal() << " : " << cosT->getVal() << " : " << evtW->getVal() << endl;
    //cout << pho1EffWeight << "\t" <<pho2EffWeight << "\t" << evtW->getVal()<<endl;

    if(evtW->getVal() >10.){
      std::cout << "Large Event Weight: " << evtW->getVal() << std::endl;
      set.Print("V");
    }
    dataMap[cat]->add(set);
  }
  cout << "Processed " << iEntry << " Entries" <<endl;

  totEB->setVal(nEB);
  totEE->setVal(nEE);

  //build a combined DataSet using the cat RooCategory
  RooArgSet setCat(set);
  setCat.add(*cat);
  RooDataSet* dataComb = new RooDataSet(tag+"_Combined","",setCat);

  float lumiRescaleFactor = (isGlobe ? 1 : lumi * 50.58/N);  //50.58 Higgs/fb

  float puWeightCorrection = (isGlobe ? 1 : (nEB+nEE)/Nentries); // globe doesn't need this correction
  
  std::cout << tag << " LUMI RESCALE: " << lumiRescaleFactor <<std::endl;

  std::map<int, RooDataSet*>::iterator dIt;
  for(dIt = dataMap.begin();dIt!=dataMap.end();dIt++){
    //loop over the barrel datasets and add the individual data to the combined dataset
    //dIt->setWeightVar(*evtW);
    TString cattag;
    cattag = Form("cat%d", dIt->first );
    std::cout << cattag <<std::endl;
    RooDataSet *tmp = new RooDataSet("DataCat_"+cattag,"",set,RooFit::Index(*cat),RooFit::Import(cattag,*(dIt->second)) );
    tmp->Print();
    Long64_t iEntry=-1;
    const RooArgSet *set;
    while( (set = tmp->get(++iEntry)) ){
      if(!isData) ((RooRealVar*)set->find("evtWeight"))->setVal( ((RooRealVar*)set->find("evtWeight"))->getVal()*lumiRescaleFactor/puWeightCorrection );
      dataComb->add(*set);
    }
    //dataComb->append(*tmp);
    //ws->import(*(dIt->second));
  }

  RooDataSet *dataComb_w = new RooDataSet(tag+"_Combined","",dataComb,setCat,0,"evtWeight");
  //import everything
  ws->import(*dataComb_w);
  ws->import(*totEB);
  ws->import(*totEE);
  ws->import(*totGen);
  ws->import(*evtW);
  cout << "Done" <<endl;

  if(efficiencyCorrection) efficiencyCorrection->Close();
  if(fileKFactor) fileKFactor->Close();
  if(fileRescaleFactor) fileRescaleFactor->Close();
  if(h) delete h;
  if(g) delete g;

  //if(f) f->Close();

}

float MakeSpinWorkspace::getEffWeight(float eta, float pt, float phi, float r9){

  if(effMaps.size()==0) return 1;

  TH3F* effMap=0;
  const char* prefix = (requireCiC ? "cic" : "pre");
  for(int i=0;i<effMaps.size(); i++){
    //std::cout << i <<std::endl;
    TObjArray *o = TString(effMaps.at(i)->GetName()).Tokenize("_");
    //the keys have names of the form <effType>_minPt_maxPt
    if( strcmp(prefix,o->At(0)->GetName()) == 0 //name match
	&& atof(o->At(1)->GetName()) <= pt  //if the min pt is below the photon pt
	&& atof(o->At(2)->GetName()) > pt){ // if the max pt is above the photon pt
      effMap = (TH3F*)effMaps.at(i);
      delete o;
      break;
    }
    delete o;
  }
  //std::cout << "found map " << effMap <<std::endl;
  float returnVal=1;
  if(effMap!=0){
    if(effMap->GetBinContent(effMap->FindFixBin(eta,phi,r9))==0){ // try to interpolate the above and below in phi
      //  std::cout << "Interpolating" <<std::endl;
      int ix,iy,iz;
      int index = effMap->FindFixBin(eta,phi,r9); // global coordinate
      effMap->GetBinXYZ(index,ix,iy,iz);
      int indexYup   = (iy == effMap->GetNbinsY() ? 1 : iy+1); //wrap around if we are at the top edge of the map
      int indexYdown = (iy == 1 ? effMap->GetNbinsY() : iy-1); //wrap around if we are at the bottom edge of the map
      
      float up = effMap->GetBinContent(ix,indexYup,iz);
      float down = effMap->GetBinContent(ix,indexYdown,iz);
      if(up==0){
	if(down!=0) returnVal= 1./down;
	else returnVal= 1;
      }
      if(down==0) returnVal= 1./up;
      returnVal= 2./(up+down); //return 1/avg
    }else{
      returnVal = 1./effMap->GetBinContent(effMap->FindFixBin(eta,phi,r9));
    }
  }
  //std::cout << "returning: " <<returnVal;
  //delete effMap;
  //delete keyList;
  if(returnVal!=returnVal){
    std::cout << "INVALID WEIGHT!!!!" <<std::endl;
    assert(false);
  }
  return returnVal;

}

void MakeSpinWorkspace::addFile(TString fName,TString l,bool is,int N,bool list){
  fileName.push_back(fName);
  label.push_back(l);
  isData.push_back(is);
  Ngen.push_back(N);
  isList.push_back(list);
  if(!is) labels->defineType(l,labels->numBins(""));
}


void MakeSpinWorkspace::MakeWorkspace(){
  //run AddToWorkspace(...) on all the input files
  std::vector<TString>::const_iterator fnIt,lIt;
  std::vector<bool>::const_iterator idIt,ilIt;
  std::vector<int>::const_iterator NgenIt;
  fnIt = fileName.begin();
  lIt  = label.begin();
  idIt = isData.begin();
  NgenIt = Ngen.begin();
  ilIt = isList.begin();

  if(!useR9 && !takeCatFromTree) getSelectionMap(selectionMap,1); // load the selection map

  for(; fnIt != fileName.end(); fnIt++, lIt++, idIt++,NgenIt++, ilIt++){
    AddToWorkspace(*fnIt,*lIt,*idIt,*NgenIt,*ilIt);
    outputFile->cd();
    ws->Write(ws->GetName(),TObject::kWriteDelete);
  }
  ws->import(*labels);

  if(mixer) mixer->mixAll();
  
  outputFile->cd();
  ws->Write(ws->GetName(),TObject::kWriteDelete);
  outputFile->Close();
}

float MakeSpinWorkspace::getEfficiency(HggOutputReader2 &h, int massPoint){
  float KFac = getKFactor(h,massPoint);
  float rescale = getRescaleFactor(h);
   if(KFac!=KFac || KFac==0){
    std::cout << "INVALID KFACTOR!!!!"<<std::endl;
    assert(false);
  }
   if(rescale!=rescale ||rescale==0){
    std::cout << "INVALID rescale!!!!"<<std::endl;
    assert(false);
  }
 return KFac*rescale;
}

float MakeSpinWorkspace::getKFactor(HggOutputReader2 &h, int massPoint){
  if(fileKFactor==0) return 1;
  TString name = Form("kfact%d_0",massPoint);

  TH1D* hist = (TH1D*)fileKFactor->Get(name);

  float kf=1;
  assert(hist!=0);
  if(hist){
    int bin = hist->FindFixBin(h.genHiggsPt);
    if(bin > hist->GetNbinsX()) bin = hist->GetNbinsX();
    if(bin <= 1) bin=2;
    kf = hist->GetBinContent(bin);
    if(kf==0) std::cout << bin << " " << h.genHiggsPt << std::endl;
  }
  delete hist;
  return kf;
}

float MakeSpinWorkspace::getRescaleFactor(HggOutputReader2 &h){
  if(fileRescaleFactor==0) return 1;

  float EFF=1;

  int cicCat = 0;
  if(requireCiC) 2*(fabs(h.PhotonPFCiC_etaSC[0])>1.48 || fabs(h.PhotonPFCiC_etaSC[1])>1.48) + (h.PhotonPFCiC_r9[0]<0.94 || h.PhotonPFCiC_r9[1]<0.94);
  else cicCat =  2*(fabs(h.Photon_etaSC[0])>1.48      || fabs(h.Photon_etaSC[1])>1.48)      + (h.Photon_r9[0]<0.94      || h.Photon_r9[1]<0.94);

  TLorentzVector p4_1; 
  TLorentzVector p4_2;

  if(requireCiC){
    p4_1.SetPtEtaPhiM(h.PhotonPFCiC_pt[0],h.PhotonPFCiC_eta[0],h.PhotonPFCiC_phi[0],0);
    p4_2.SetPtEtaPhiM(h.PhotonPFCiC_pt[1],h.PhotonPFCiC_eta[1],h.PhotonPFCiC_phi[1],0);

  }else{
    p4_1.SetPtEtaPhiM(h.Photon_pt[0],h.Photon_eta[0],h.Photon_phi[0],0);
    p4_2.SetPtEtaPhiM(h.Photon_pt[1],h.Photon_eta[1],h.Photon_phi[1],0);
  }

  float pt_sys = (p4_1+p4_1).Pt();
  
  //L1-HLT factor
  TGraphAsymmErrors* e = (TGraphAsymmErrors*)fileRescaleFactor->Get( Form("effL1HLT_cat%d",cicCat) );
  EFF*=getEffFromTGraph(e,pt_sys);
  delete e;
  //vertex efficiencies
  e = (TGraphAsymmErrors*)fileRescaleFactor->Get( Form("ratioVertex_cat%d_pass",cicCat) );
  EFF*=getEffFromTGraph(e,pt_sys);
  delete e;
  e = (TGraphAsymmErrors*)fileRescaleFactor->Get( Form("ratioVertex_cat%d_fail",cicCat) );
  EFF*=getEffFromTGraph(e,pt_sys);
  delete e;
  //per-photon efficiencies
  TString labels[4] = {"EBHighR9","EBLowR9","EEHighR9","EELowR9"};
  for(int i=0;i<2;i++){
    int index = 2*(fabs(h.Photon_etaSC[i])>1.48)+(h.Photon_r9[i] > 0.94);
    e = (TGraphAsymmErrors*)fileRescaleFactor->Get( TString("ratioTP_")+labels[index] );
    if(requireCiC) EFF*=getEffFromTGraph(e,h.PhotonPFCiC_pt[i]);
    else EFF*=getEffFromTGraph(e,h.Photon_pt[i]);
    delete e;
    e = (TGraphAsymmErrors*)fileRescaleFactor->Get( TString("ratioR9_")+labels[index] );
    if(requireCiC) EFF*=getEffFromTGraph(e,h.PhotonPFCiC_pt[i]);
    else EFF*=getEffFromTGraph(e,h.Photon_pt[i]);
    delete e;
  }
  return EFF;
}


float MakeSpinWorkspace::getEffFromTGraph(TGraphAsymmErrors* e,float pt){
  int nBins = e-> GetN();
  int selectedBin=0;
  for(selectedBin=0;selectedBin<nBins;selectedBin++){
    Double_t x,y;
    e->GetPoint(selectedBin,x,y);
    if(pt < x) break; // found the first bin above the pT
  }
  if(selectedBin==nBins) selectedBin--;                       

  if(selectedBin==0){
    double x,y;
    e->GetPoint(selectedBin,x,y);
    return y;
  }else{
    double xLow,yLow;
    e->GetPoint(selectedBin-1,xLow,yLow);
    double xHigh,yHigh;
    e->GetPoint(selectedBin,xHigh,yHigh);
    float fracLow = (pt-xLow)/(xHigh-xLow); // the weight to give to the lower bin
    float weight = yLow*fracLow+yHigh*(1-fracLow);
    return weight;
  }
  return 1;
}

float MakeSpinWorkspace::calculateCosThetaCS(HggOutputReader2 *h){
  TLorentzVector p4_1; p4_1.SetPtEtaPhiM(h->Photon_pt[0],h->Photon_eta[0],h->Photon_phi[0],0);
  TLorentzVector p4_2; p4_2.SetPtEtaPhiM(h->Photon_pt[1],h->Photon_eta[1],h->Photon_phi[1],0);

  TLorentzVector f = p4_1+p4_2;
  TVector3 fBoost = -1*f.BoostVector();

  TLorentzVector b1; b1.SetPxPyPzE(0,0,4000,4000);
  TLorentzVector b2; b2.SetPxPyPzE(0,0,-4000,4000);
  //b1.Print();
  b1.Boost(fBoost);
  b2.Boost(fBoost);
  p4_1.Boost(fBoost);
  p4_2.Boost(fBoost);
  //b1.Print();
  //return TMath::Sqrt(1 - f.Pt()*f.Pt()/( f.Pt()*f.Pt() + f.M2() ));

  TVector3 g1_unit = p4_1.Vect().Unit();
  TVector3 b1_unit = b1.Vect().Unit();
  TVector3 b2_unit = b2.Vect().Unit();
  TVector3 direction_cs = (b1_unit - b2_unit).Unit();
  
  return TMath::Abs(TMath::Cos(direction_cs.Angle(g1_unit)));

  //return TMath::Cos((TMath::Pi()- b1.Angle(b2.Vect()))/2);
}


TChain* MakeSpinWorkspace::getChainFromList(TString inputFileList, TString treeName){
  TChain *chain = new TChain(treeName,treeName);

  std::ifstream listFile(inputFileList.Data());

  std::string line;
  while(std::getline(listFile,line)){
    if(line.size() < 2) continue;
    chain->AddFile(line.c_str());
  }
  return chain;
}

bool MakeSpinWorkspace::passCiCIso(HggOutputReader2 &h, int i){

  int index = (fabs(h.Photon_etaSC[i]) > 1.48)*2+(h.Photon_r9[i]<0.94);

  if(h.Photon_dr03PFChargedIso[i] > chargedIso[index]) return false;
  if(h.Photon_isosumGood[i] > goodIsoSum[index]) return false;
  if(h.Photon_isosumBad[i] > badIsoSum[index]) return false;

  return true;
}

void MakeSpinWorkspace::setupBranches(HggOutputReader2 &h){
  h.fChain->SetBranchStatus("*",0);

  h.fChain->SetBranchStatus("nPhoton",1);
  h.fChain->SetBranchStatus("Photon.*",1);
  h.fChain->SetBranchStatus("PhotonPFCiC.*",1);
  h.fChain->SetBranchStatus("mPair",1);
  h.fChain->SetBranchStatus("mPairPFCiC",1);
  h.fChain->SetBranchStatus("mPairNoCorr",1);
  h.fChain->SetBranchStatus("mPairNoCorrPFCiC",1);
  h.fChain->SetBranchStatus("cosThetaLead",1);
  h.fChain->SetBranchStatus("cosThetaLeadPFCiC",1);
  h.fChain->SetBranchStatus("evtWeight",1);
  h.fChain->SetBranchStatus("diPhotonMVA",1);
  h.fChain->SetBranchStatus("runNumber",1);
  h.fChain->SetBranchStatus("genHiggsPt",1);
}
