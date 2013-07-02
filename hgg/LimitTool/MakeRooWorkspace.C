//for RooFit:
#include "RooRealVar.h"               
#include "RooArgSet.h"                
#include "RooArgList.h"               
#include "RooDataSet.h"               
#include "RooExponential.h"           
#include "RooLandau.h"                
#include "RooPlot.h"                  
#include "RooFit.h"                   
#include "RooAddPdf.h"                
#include "RooGaussian.h"              
#include "RooCBShape.h"               
#include "RooFFTConvPdf.h"            
#include "RooDataHist.h"              
#include "RooHistPdf.h"               
#include "RooHistFunc.h"              
#include "RooMoment.h"                
#include "RooFitResult.h"             
#include "RooExtendPdf.h"             
#include "RooGenericPdf.h"            
#include "RooBreitWigner.h"           
#include "RooBifurGauss.h"            
#include "RooProdPdf.h"               
#include "RooCategory.h"              
#include "RooSimultaneous.h"          
#include "RooWorkspace.h"             
#include "RooConstVar.h"              
#include "TEfficiency.h"              
#include "RooConstVar.h"              
#include "RooKeysPdf.h"               
#include "RooBernstein.h"             
#include "RooPolynomial.h"            
#include "RooIntegralMorph.h"         
#include "RooNLLVar.h"                
#include "RooAbsArg.h"

#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <TList.h>
#include <TObjArray.h>

#include "defineCategories.C"

using namespace std;
using namespace RooFit;               

#define debug 0
#include "branchingRatio.cc"
#define nCat_ 6
struct WorkspaceVars{
  RooWorkspace * outputWorkspace;
  std::map<std::string, RooRealVar*>  *m_real_var_;
  std::map<std::string, double> *m_var_min_;
  std::map<std::string, double> *m_var_max_;  
  std::map<std::string,RooDataSet*> *data_;
  float massMin_,massMax_;
};

void createDataSet(WorkspaceVars*,string,string);
void addDataPoint(WorkspaceVars*,string, float, float);
double getHistVal(TH1D* hist, float x);
vector<float> getKFac(std::vector<TH1D*>* kFacs,float higgsPt,int shift);

int getPhotonCat(float eta,float r9){
  return (r9<0.94)+2*(fabs(eta)>1.48);
}

typedef std::pair<string,std::pair<int,int> > indexType;
typedef std::map<std::string,std::vector<TGraphAsymmErrors*> > EffMap;
EffMap* loadEffCor(TFile* f);
vector<float> effCorrMC(float pt, TString type, int cat, EffMap* map,int shift=0);


TObjArray* MakeRooWorkspace(RooWorkspace *workspace,TString fileName, bool isData,float lumi,int massPoint,TString process, TString kFactorFile="",TString pileupReWeightFName="",TString EffCorrectionFName="",TString selection="MVA"){
  cout << "Making workspace for mass point: " << massPoint <<  "  Process: " << process << endl;
  if(debug) cout << "MakeRooWorkspace" <<endl;
  TChain *chain = new TChain("HggOutput");
  chain->AddFile(fileName);


  
  //roohist_sig_ggh_mass_m$MASS_$CHANNEL
  TObjArray *outArray = new TObjArray();
  const int    nSysts  = 7;
  const string systs[] = {"E_res","E_scale","idEff","kFactor","r9Eff","triggerEff","vtxEff"};

  std::map< indexType, int> nameIndexMap;

  const int nSigma=1;
  int ind=0;

  std::map<std::string,std::vector<float> > efficiencies; // map between the efficiency name and the vector of efficiencies

  if(!isData){
    for(int iCat=0;iCat<nCat_;iCat++){
      outArray->Add(new TH1F(Form("th1f_sig_%s_mass_m%d_cat%d",process.Data(),massPoint,iCat),"",160,100,180));
      nameIndexMap[indexType(Form("cat%d",iCat),std::pair<int,int>(0,iCat))] = ind++;
      for(int iSyst=0;iSyst<nSysts;iSyst++){
	for(int iSig= -nSigma; iSig<=nSigma; iSig++){
	  if(iSig==0) continue;
	  outArray->Add(new TH1F(Form("th1f_sig_%s_mass_m%d_cat%d_%s%s%02d_sigma",process.Data(),massPoint,iCat,systs[iSyst].c_str(),(iSig>0?"Up":"Down"),abs(iSig)),"",160,100,180));
	  if(debug) cout << Form("th1f_sig_%s_mass_m%d_cat%d_%s%s%02d_sigma",process.Data(),massPoint,iCat,systs[iSyst].c_str(),(iSig>0?"Up":"Down"),abs(iSig)) << endl;
	  nameIndexMap[indexType(systs[iSyst],std::pair<int,int>(iSig,iCat))] = ind++;
	}
      }
    
    }
  }
  WorkspaceVars vars;
  vars.m_real_var_ = new std::map<std::string, RooRealVar*>;
  vars.m_var_min_ = new std::map<std::string, double>;
  vars.m_var_max_ = new std::map<std::string, double>;  
  vars.data_ = new std::map<std::string,RooDataSet*>;
  vars.massMin_=100;
  vars.massMax_=180;

  TFile * pileupReWeightFile = new TFile(pileupReWeightFName);
  TH1F* pileupReWeight = (TH1F*)pileupReWeightFile->Get("pileupReWeight");

  TFile *EffCorrFile = new TFile(EffCorrectionFName);
  EffMap* EfficiencyMap = loadEffCor(EffCorrFile);
  if(debug) cout << "Loaded All Maps" << endl;
  float mPair;
  float MVA;
  float genHiggsPt;
  float genHiggsVx;
  float genHiggsVy;
  float genHiggsVz;
  float nPU;
  float Vx;
  float Vy;
  float Vz;
  float ptPho[2];
  float etaPho[2];
  float phiPho[2];
  float energyPho[2];
  float r9Pho[2];
  float Mjj;
  vector<float> *mPairScale = 0;
  vector<float> *mPairSmear = 0;
  chain->SetBranchAddress("mPair",&mPair);
  chain->SetBranchAddress("diPhotonMVA",&MVA);
  chain->SetBranchAddress("genHiggsPt",&genHiggsPt);
  chain->SetBranchAddress("genHiggsVx",&genHiggsVx);
  chain->SetBranchAddress("genHiggsVy",&genHiggsVy);
  chain->SetBranchAddress("genHiggsVz",&genHiggsVz);
  chain->SetBranchAddress("nPU",&nPU);
  chain->SetBranchAddress("diPhotonVtxX",&Vx);
  chain->SetBranchAddress("diPhotonVtxY",&Vy);
  chain->SetBranchAddress("diPhotonVtxZ",&Vz);
  chain->SetBranchAddress("Photon.pt",ptPho);
  chain->SetBranchAddress("Photon.eta",etaPho);
  chain->SetBranchAddress("Photon.phi",phiPho);
  chain->SetBranchAddress("Photon.E",energyPho);
  chain->SetBranchAddress("Photon.r9",r9Pho);
  chain->SetBranchAddress("Mjj",&Mjj);

  chain->SetBranchAddress("mPairScale",&mPairScale); //for E_scale systematic
  chain->SetBranchAddress("mPairSmear",&mPairSmear); //for E_res   systematic

  CatMap categories = getCategoryCuts(chain);

  
  branchingRatio XSECS;

  std::vector<TH1D*> kFacs;
  bool doKFactor = false;
   if(kFactorFile!=""){
    TFile *file = new TFile(kFactorFile);
    /*
    kFacs.push_back((TH1D*)file->Get(Form("kfact%d_0",int(massPoint)) )); //kFac
    kFacs.push_back((TH1D*)file->Get(Form("kfact%d_1",int(massPoint)))); //UP
    kFacs.push_back((TH1D*)file->Get(Form("kfact%d_6",int(massPoint)))); //DOWN
    */
    kFacs.push_back((TH1D*)file->Get(Form("kfact%d_0",120) )); //kFac
    kFacs.push_back((TH1D*)file->Get(Form("kfact%d_1",120))); //UP
    kFacs.push_back((TH1D*)file->Get(Form("kfact%d_6",120))); //DOWN

    doKFactor = true;
  }

   //RooWorkspace outputWorkspace;
  int nMassBins_;
  //int nCat_;
  bool isSigMC = true;


  //setup the RooWorkspace
  //outputWorkspace
  std::map<std::string, RooRealVar*> test;
  if(debug) cout << "Creating Workspace" <<endl;
  RooRealVar *tmp = new RooRealVar("mass","mass",0.,vars.massMin_,vars.massMax_);
  if(debug) cout << "Inserting mass var" << endl;
  vars.m_real_var_->insert(std::pair<string,RooRealVar*>(string("mass"),tmp));
  tmp = new RooRealVar("weight_mass","weight_mass",0.,0,1e6);
  if(debug) cout << "Inserting mass weigth var" << endl;
  vars.m_real_var_->insert(std::pair<string,RooRealVar*>(string("weight_mass"),tmp));
  createDataSet(&vars,"mass","data_mass");
  createDataSet(&vars,"mass","bkg_mass");
  if(isSigMC){
    string dname = string(Form("sig_%s_mass_m%d",process.Data(),massPoint));
    createDataSet(&vars,"mass",dname);
    dname = string(Form("sig_%s_mass_rv_m%d",process.Data(),massPoint));
    createDataSet(&vars,"mass",dname);
    dname = string(Form("sig_%s_mass_wv_m%d",process.Data(),massPoint));
    createDataSet(&vars,"mass",dname);
  }
  if(debug) cout << "Done" << endl << "Looping" << endl;
  cout << "Getting Entries" << endl;
  long nEntries = chain->GetEntries();
  long entry=-1;
  Int_t TreeNum = -99;
  cout << "Looping, Entries: "<< nEntries << endl;
  while(chain->GetEntry(++entry)){
    if(entry%500==0) cout << "Processing Entry: " << entry << "\r" << flush;
    if(chain->GetTreeNumber() != TreeNum){
      updateFormulas(&categories);
      TreeNum = chain->GetTreeNumber();
    }
    int cat = getCat(&categories,selection,nCat_);
    if(debug) cout << "cat: " << cat << endl;
    if(cat<0) continue;
    if(isData) addDataPoint(&vars,Form("data_mass_cat%d",cat),mPair,1);
    else{ // not data!!
      TVector3 genVtx(genHiggsVx,genHiggsVy,genHiggsVz);
      TVector3 recoVtx(Vx,Vy,Vz);
      TLorentzVector p1_p4; p1_p4.SetPtEtaPhiM(ptPho[0],etaPho[0],phiPho[0],energyPho[0]);
      TLorentzVector p2_p4; p2_p4.SetPtEtaPhiM(ptPho[1],etaPho[1],phiPho[1],energyPho[1]);
      TLorentzVector sys_p4 = p1_p4+p2_p4;
      bool CorrectVtx = (genVtx-recoVtx).Mag() < 1.;
      if(debug) cout << "CorrectVtx: " << CorrectVtx << endl;
      efficiencies["kFactor"] = getKFac(&kFacs,genHiggsPt,nSigma);
      if(debug) cout << "kFactor: " << efficiencies["kFactor"].at(nSigma+1) << endl;
      float XSweight = XSECS.BranchingRatios[massPoint]*(XSECS.Proc8TeV[process])[massPoint]*lumi*1000./nEntries;
      if(debug) cout << "XSweight: " << XSweight << endl; 
      float puwt=pileupReWeight->GetBinContent( pileupReWeight->FindFixBin(nPU) );
      if(debug) cout << "puwt: " << puwt << endl;

      float photonEff=1;
      std::vector<float> tpPho1 = effCorrMC(p1_p4.Pt(),"ratioTP",getPhotonCat(etaPho[0],r9Pho[0]),EfficiencyMap,nSigma);
      std::vector<float> tpPho2 = effCorrMC(p2_p4.Pt(),"ratioTP",getPhotonCat(etaPho[1],r9Pho[1]),EfficiencyMap,nSigma);
      std::vector<float> r9effPho1 = effCorrMC(p1_p4.Pt(),"ratioR9",getPhotonCat(etaPho[0],r9Pho[0]),EfficiencyMap,nSigma);
      std::vector<float> r9effPho2 = effCorrMC(p2_p4.Pt(),"ratioR9",getPhotonCat(etaPho[1],r9Pho[1]),EfficiencyMap,nSigma);
      std::vector<float> idEff, r9Eff;
      for(int i=0;i<tpPho1.size();i++){
	idEff.push_back(tpPho1.at(i)*tpPho2.at(i));
	r9Eff.push_back(r9effPho1.at(i)*r9effPho2.at(i));
      }
      efficiencies["idEff"] = idEff;
      efficiencies["r9Eff"] = r9Eff;

      float ptPair = sys_p4.Pt();
      if(debug) cout << "ptPair: " << ptPair << endl;
      float diPhotonEff=1;
      efficiencies["triggerEff"] = effCorrMC(ptPair,"effL1HLT",cat,EfficiencyMap,nSigma);
      if(CorrectVtx) efficiencies["vtxEff"] = effCorrMC(ptPair,"ratioVertex_pass",cat,EfficiencyMap,nSigma);
      else           efficiencies["vtxEff"] = effCorrMC(ptPair,"ratioVertex_fail",cat,EfficiencyMap,nSigma);

      int ci = nSigma+1; // this is the index of the central value
      float evtWeight = efficiencies["kFactor"].at(ci)*efficiencies["idEff"].at(ci)*efficiencies["r9Eff"].at(ci);
      evtWeight*=efficiencies["triggerEff"].at(ci)*efficiencies["vtxEff"].at(ci)*XSweight*puwt;

      addDataPoint(&vars,Form("sig_%s_mass_m%d_cat%d",process.Data(),massPoint,cat),mPair,evtWeight);
      if(CorrectVtx) addDataPoint(&vars,Form("sig_%s_mass_rv_m%d_cat%d",process.Data(),massPoint,cat),mPair,evtWeight);
      else addDataPoint(&vars,Form("sig_%s_mass_wv_m%d_cat%d",process.Data(),massPoint,cat),mPair,evtWeight);

      ((TH1F*)outArray->At(nameIndexMap[indexType(string(Form("cat%d",cat)),std::pair<int,int>(0.,cat))]))->Fill(mPair,evtWeight);
      if(debug) cout << "Filling Systematics:" << endl;
       //SYSTEMATICS!!!!
      for(int iSyst=0; iSyst < nSysts; iSyst++){
	if(debug) cout << "This Syst: " << systs[iSyst] << endl;
	for(int iSigma = -nSigma;iSigma<=nSigma;iSigma++){
	  if(iSigma==0) continue;
	  if(iSyst<2){ //energy scale/smear
	    float thisEvtWeight = evtWeight;  // use the default weighting
	    float thisMPair = mPair;
	    if(iSyst==0) thisMPair = mPairSmear->at(iSigma+nSigma);
	    else         thisMPair = mPairScale->at(iSigma+nSigma);
	    if(debug) cout << "This EvtWeight: " << thisEvtWeight << endl 
			   << "thisMass: " << thisMPair << "   default: " << mPair << endl;
	    ((TH1F*)outArray->At(nameIndexMap[indexType(systs[iSyst],std::pair<int,int>(iSigma,cat))]))->Fill(thisMPair,thisEvtWeight);	    
	  }else{ //these are event weight scalings
	    if(debug){
	      cout << "this Sigma: " << efficiencies[systs[iSyst]].at(iSigma+nSigma) << endl;
	      cout << "def  Sigma: " << efficiencies[systs[iSyst]].at(ci)     << endl;
	    }
	    float thisEvtWeight = evtWeight*efficiencies[systs[iSyst]].at(iSigma+nSigma)/efficiencies[systs[iSyst]].at(ci);
	    ((TH1F*)outArray->At(nameIndexMap[indexType(systs[iSyst],std::pair<int,int>(iSigma,cat))]))->Fill(mPair,thisEvtWeight);
	  }
	}
      }//for(iSyst...
      if(mPair > 50 && debug) exit(0);
      
      
    } //!isData
  }
  std::map<std::string,RooDataSet*>::iterator it_data;
  for(it_data = vars.data_->begin();it_data!=vars.data_->end();it_data++) {
    workspace->import(*(it_data->second));
  } 
  //outArray->Write();
  EffCorrFile->Close();
  //TFile *f = new TFile(outputName,"RECREATE");
  //outputWorkspace.Write();
  //f->Close();
  delete vars.m_real_var_;
  delete vars.m_var_min_;
  delete vars.m_var_max_;
  delete vars.data_;
  cout << endl;
  return outArray;
  }

void createDataSet(WorkspaceVars* vars,string name,string data_name_base){
  if(debug) cout << "createDataSet" << endl;
  for(int i=0;i<nCat_;i++){
    if(debug) cout << i << endl;
    //createDataSet("mass","mass_cat%d",nbins,x1,x2)
    //std::map<std::string,RooRealVar>::const_iterator test=m_real_var_.find(name);
    string data_name = string(Form("%s_cat%d",data_name_base.c_str(),i));
    if(debug) cout << "data_name: " << data_name << endl;
    vars->m_var_min_->insert(pair<std::string, double >(data_name,vars->massMin_)); 
    vars->m_var_max_->insert(pair<std::string, double >(data_name,vars->massMax_));
    RooRealVar *mass = (*(vars->m_real_var_))["mass"];
    RooRealVar *mass_w = (*(vars->m_real_var_))["weight_mass"];
    if(debug) cout << "mass " << mass->getVal() << endl;
    if(debug) cout << "mass_w " << mass_w->getVal() << endl;
    RooArgSet *args = new RooArgSet(*mass,*mass_w);
    //RooRealVar test("mass","mass",0,-10,10);
    //RooArgSet args(test);
    if(debug) cout << "making data set" << endl;
    RooDataSet *data_tmp = new RooDataSet(data_name.c_str(),data_name.c_str(),*args,WeightVar("weight_mass"));
    
    //RooDataSet *data_tmp = new RooDataSet(data_name.c_str(),data_name.c_str(),args,"mass");
    if(debug) cout << "insert" << endl;
    vars->data_->insert(std::pair<std::string,RooDataSet*>(data_name,data_tmp));
  }
}

void addDataPoint(WorkspaceVars* vars,string name, float mass, float weight){
  if(debug) cout << "addDataPoint m=" << mass << "  name=" << name << endl;
  if(mass < vars->massMin_ || mass > vars->massMax_) return;
  if(debug) cout << "Setting Var" << endl;
  *(*(vars->m_real_var_))["mass"] = mass;
  //*(*(vars->m_real_var_))["weight_mass"] = weight;
  if(debug) cout << "Adding" << endl;
  (*(vars->data_))[name]->add(RooArgSet( *(*(vars->m_real_var_))["mass"]),weight);
  //string name = string(Form("data_mass_cat%d",cat));
  //data_[name].add( RooArgSet( RooRealVar(mass) ), weight);
}

vector<float> getKFac(std::vector<TH1D*>* kFacs,float higgsPt,int shiftRange){


  float mean = getHistVal(kFacs->at(0),higgsPt);
  float up = getHistVal(kFacs->at(1),higgsPt);
  float down = getHistVal(kFacs->at(2),higgsPt);
  std::vector<float> out;
  for(int iShift=-shiftRange; iShift <= shiftRange; iShift++){
    if(iShift<0) out.push_back(mean + (down-mean)*fabs(iShift));
    if(iShift>0) out.push_back(mean + (up  -mean)*fabs(iShift));
    else         out.push_back(mean);
  }
  return out;
}

double getHistVal(TH1D* hist, float x){ // returns the value corresponding to a given x (returns the top of bottom bin rather than over/underflow
  if(hist==0) cout << "ERROR: TRYING TO GET VALUE FROM AN INVALID HISTOGRAM" << endl;
  if(x < hist->GetXaxis()->GetXmin()) return hist->GetBinContent(1);
  int nBins = hist->GetNbinsX();
  if(x > hist->GetXaxis()->GetXmax()) return hist->GetBinContent(nBins);
  return hist->GetBinContent( hist->FindFixBin(x) ); 
}


EffMap* loadEffCor(TFile* f){
  EffMap* out = new EffMap;
  const int nNames=5;
  const string labels[] = {"ratioTP","ratioR9","ratioVertex_pass","ratioVertex_fail","effL1HLT"};
  char *names[] = {
    "ratioTP_%s",
    "ratioR9_%s",
    "ratioVertex_cat%d_pass",
    "ratioVertex_cat%d_fail",
    "effL1HLT_cat%d"
  };
 
  //const int nCat = 4;
  const string catNames[] = {"EBHighR9","EBLowR9","EEHighR9","EELowR9"};
  //std::map<std::string,std::vector<TGraphAsymmErrors*> >
  
  for(int i=0;i<nNames;i++){
    std::vector<TGraphAsymmErrors*> tmp;
    for(int j=0;j<nCat_;j++){
      if(i>=2)  tmp.push_back( (TGraphAsymmErrors*)f->Get( Form(names[i],j) ) );
      else     tmp.push_back( (TGraphAsymmErrors*)f->Get( Form(names[i],catNames[j].c_str()) ) );
    }
    (*out)[labels[i]] = tmp;
  }
  return out;
}

vector<float> effCorrMC(float pt, TString type, int cat, EffMap* map,int shiftRange){
  if(debug) cout << "effCorrMC: pt=" << pt << " type=" << type.Data() << " cat=" << cat << endl;
  TGraphAsymmErrors* e = (map->find(type.Data()))->second.at(cat);

  int nBins = e-> GetN();
  int selectedBin=0;
  for(selectedBin=0;selectedBin<nBins;selectedBin++){
    Double_t x,y;
    e->GetPoint(selectedBin,x,y);
    if(pt < x) break; // found the first bin above the pT
  }
    if(selectedBin==nBins) selectedBin--;
  if(debug) cout << "\tSelectedBin: " << selectedBin << endl;
  vector<float> out;
  for(int iShift = -shiftRange; iShift<=shiftRange; iShift++){
    if(selectedBin==0 || selectedBin==nBins-1){ // its below the first or above the last bin
      double x,y;
      e->GetPoint(selectedBin,x,y);
      if(debug) cout << "\ty: " << y << endl;
      if(iShift==0) out.push_back(y);
      else if(iShift<0) out.push_back(y+iShift*e->GetErrorYlow(selectedBin-1));
      else out.push_back(y+iShift*e->GetErrorYhigh(selectedBin-1));
    }
  
    double xLow,yLow;
    e->GetPoint(selectedBin-1,xLow,yLow);
    double xHigh,yHigh;
    e->GetPoint(selectedBin,xHigh,yHigh);
    if(debug) cout << "\tyLow: " << yLow << "  yHigh: " << yHigh << endl;
    float fracLow = (pt-xLow)/(xHigh-xLow); // the weight to give to the lower bin
    if(debug) cout << "\tfracLow: " << fracLow << endl;
    float weight = yLow*fracLow+yHigh*(1-fracLow);
    float error = e->GetErrorYlow(selectedBin)*fracLow+e->GetErrorYhigh(selectedBin+1)*(1-fracLow);
    if(debug) cout << "weight: " << weight << endl;
    out.push_back(weight + iShift*error);
  }
  return out;
}



void MakeAllWorkspaces(string fileList, float lumi,string outputFile,
		       TString kFactorFileName, TString pileupReWeightDirName, 
		       TString EfficiencyCorrectionFileName,TString SelectionType){

  RooWorkspace *outputWorkspace = new RooWorkspace();
  outputWorkspace->SetName("cms_hgg_workspace");
  
  ifstream fileListStream(fileList.c_str(),fstream::in);

  string file;
  string outputName;
  string puName;
  string massString;
  string process;
  const int nProc = 4;
  const string processes[nProc] = {"GluGluToH","WH_ZH","ttH","VBF_H"};
  const string procNames[nProc] = {"ggh","wzh","tth","vbf"};
  TFile *f = new TFile(outputFile.c_str(),"RECREATE");

  while(fileListStream.good()){
    getline(fileListStream,file);
    if(file.empty()) continue;
    if(file.find("/castor/cern.ch")!=string::npos)
      file.insert(0,"rfio://");

    outputName = file.substr(file.find_last_of('/')+1);
    outputName.insert(outputName.find(".root"),"_WorkSpace");
    outputName.insert(0,"workspaces/");
    cout << outputName << endl;
    puName = file.substr(file.find_last_of('/')+1);
    puName.insert(puName.find(".root"),"_puReWeight");
    puName.insert(0,"/");
    puName.insert(0,pileupReWeightDirName.Data());
    cout << puName <<endl;
    //find mass point
    massString = file.substr(file.find_last_of("M")+2,3);
    cout << file << "  " << outputName << "  " << puName << "  " << massString << endl;

    for(int i=0;i<nProc;i++){
      if(file.find(processes[i])!=string::npos){
	process = procNames[i];
	break;
      }
    }

    TObjArray* hists = MakeRooWorkspace(outputWorkspace,TString(file), 0, lumi,atoi(massString.c_str()),
					TString(process),
					kFactorFileName, puName,
					EfficiencyCorrectionFileName,SelectionType);
    f->cd();
    if(hists) hists->Write();
    //for(int i=0;i<hists->GetEntries();i++) hists->At(i)->Write();
    
  }
  fileListStream.close();

  outputWorkspace->Write();
  f->Close();
}
