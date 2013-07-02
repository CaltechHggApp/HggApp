#include "MakeSpinWorkspace.C"
#include "MakeSpinPlots.C"

void MakeDataSet(TString inputFile,RooWorkspace* ws,TString tag, bool isData,bool requireCiC=false,TString mcRescaleFile="",int runMin=0, int runMax=999999){
  TFile *f = TFile::Open(inputFile);
  TTree *tree = (TTree*)f->Get("HggOutput");
  
  HggOutputReader2 h(tree);
  
  float **offset=0,**scale=0;
  if(mcRescaleFile!=""){
    TFile mcRescale(mcRescaleFile);
    getSigEoEScales(&mcRescale,offset,scale);
  }
  // fit variables
  RooRealVar* mass   = new RooRealVar("mass",  "Mass [GeV]", 100., 180.);
  RooRealVar* cosT   = new RooRealVar("cosT",  "cos(theta)", -1, 1);
  RooRealVar* sige1  = new RooRealVar("sigEoE1","#sigma_{E}/E Lead Photon",0,0.1);
  RooRealVar* sige2  = new RooRealVar("sigEoE2","#sigma_{E}/E SubLead Photon",0,0.1);
  RooRealVar* evtW   = new RooRealVar("evtWeight","Event Weight",1,0,1e6);

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
  

  RooArgSet set;
  set.add(*mass);
  set.add(*cosT);
  set.add(*sige1);
  set.add(*sige2);
  set.add(*evtW);
  set.add(*eta1); set.add(*eta2);
  set.add(*etaSC1); set.add(*etaSC2);
  set.add(*phi1); set.add(*phi2);
  set.add(*pt1); set.add(*pt2);
  set.add(*r91); set.add(*r92);
  set.add(*idMVA1); set.add(*idMVA2);
  set.add(*diPhotonMVA);

  RooRealVar *totEB = new RooRealVar(Form("%s_EB_totalEvents",tag.Data()),"",0,0,1e9);
  RooRealVar *totEE = new RooRealVar(Form("%s_EE_totalEvents",tag.Data()),"",0,0,1e9);

  std::map<std::pair<int,int>, RooDataSet*> dataMapEB, dataMapEE;
  std::map<std::pair<int,int>, RooDataSet*> *datamap;

  RooDataSet* dataEB = new RooDataSet(Form("%s_EB",tag.Data()),"",set,"evtWeight");
  RooDataSet* dataEE = new RooDataSet(Form("%s_EE",tag.Data()),"",set,"evtWeight");

  Long64_t iEntry=-1;
  cout << "Making DataSet" << endl;
  Long64_t nEB=0,nEE=0;
  while(h.GetEntry(++iEntry)){
    if( !(iEntry%10000) ) cout << "Processing " << iEntry <<endl;
    int maxI = (h.Photon_pt[1] > h.Photon_pt[0] ? 1:0);
    int minI = (h.Photon_pt[1] > h.Photon_pt[0] ? 0:1);
    if(fabs(h.Photon_etaSC[1]) < 1.48 && fabs(h.Photon_etaSC[0]) < 1.48) nEB++;
    else nEE++;

    if(!getBaselineSelection(&h,maxI,minI)) continue;
    if(isData && (h.runNumber < runMin || h.runNumber > runMax) ) continue;

    if(requireCiC){
      if(h.Photon_passPFCiC[1]==false || h.Photon_passPFCiC[0]==false) continue;
    }

    float se1 = h.Photon_EError[maxI]/h.Photon_E[maxI];
    float se2 = h.Photon_EError[minI]/h.Photon_E[minI];
    if(!isData){
      
      se1 = scaleSEoE(se1,h.Photon_etaSC[maxI],h.Photon_pt[maxI],offset,scale);
      se2 = scaleSEoE(se2,h.Photon_etaSC[minI],h.Photon_pt[minI],offset,scale);
    }
    mass->setVal(h.mPair);
    cosT->setVal(h.cosThetaLead);
    sige1->setVal(se1);
    sige2->setVal(se2);
    evtW->setVal(1);

    eta1->setVal(h.Photon_eta[maxI]);
    eta2->setVal(h.Photon_eta[minI]);
    
    etaSC1->setVal(h.Photon_etaSC[maxI]);
    etaSC2->setVal(h.Photon_etaSC[minI]);
    
    phi1->setVal(h.Photon_phi[maxI]);
    phi2->setVal(h.Photon_phi[minI]);
    
    pt1->setVal(h.Photon_pt[maxI]);
    pt2->setVal(h.Photon_pt[minI]);
    
    r91->setVal(h.Photon_r9[maxI]);
    r92->setVal(h.Photon_r9[minI]);
    
    idMVA1->setVal(h.Photon_idMVA[maxI]);
    idMVA2->setVal(h.Photon_idMVA[minI]);
    
    diPhotonMVA->setVal(h.diPhotonMVA);

    if( fabs(h.Photon_etaSC[1]) < 1.48 && fabs(h.Photon_etaSC[0]) < 1.48 ) dataEB->add(set);
    else dataEE->add(set);
  }
  cout << "Processed " << iEntry << " Entries" <<endl;

  totEB->setVal(nEB);
  totEE->setVal(nEE);

  std::map<std::pair<int,int>, RooDataSet*>::iterator dIt;
  for(dIt = dataMapEB.begin();dIt!=dataMapEB.end();dIt++) ws->import(*(dIt->second));
  for(dIt = dataMapEE.begin();dIt!=dataMapEE.end();dIt++) ws->import(*(dIt->second));
  ws->import(*dataEB);
  ws->import(*dataEE);

  ws->import(*totEB);
  ws->import(*totEE);
  cout << "Done" <<endl;
  if(offset) delete offset;
  if(scale)  delete scale;
}


void MakeOptimizationWorkspace(TString inputDataFile, TString inputHggFile,TString outputFile,bool requireCiC=false){
  TFile *out = new TFile(outputFile,"RECREATE");
  RooWorkspace * ws = new RooWorkspace();
  ws->SetName("cms_hgg_spin_total_workspace");  

  MakeDataSet(inputDataFile,ws,"Data",true,requireCiC,"");
  MakeDataSet(inputHggFile,ws,"Hgg125",true,requireCiC,"");
    
  out->cd();
  ws->Write();
  out->Close();
}

void MakeCutWorkspaces(RooWorkspace* inputWs,RooWorkspace* outputWs,TString tag,float offsetEB,float offsetEE){
  // parameters: [0]+exp([1]+[2]*pt;
  //

  TH2F* map = getSelectionMap5(true);
  //offset the map
  for(int i=0;i<map->GetNbinsY();i++){
    map->SetBinContent(1,i,map->GetBinContent(1,i)+offsetEB);
    map->SetBinContent(2,i,map->GetBinContent(2,i)+offsetEB);
    map->SetBinContent(3,i,map->GetBinContent(3,i)+offsetEE);
    map->SetBinContent(4,i,map->GetBinContent(4,i)+offsetEE);
  }

  RooDataSet* EB = (RooDataSet*)inputWs->data(Form("%s_EB",tag.Data()));
  RooDataSet* EE = (RooDataSet*)inputWs->data(Form("%s_EE",tag.Data()));

  RooRealVar *totEB = new RooRealVar(Form("%s_EB_totalEvents",tag.Data()),"",EB->sumEntries(),0,1e9);
  RooRealVar *totEE = new RooRealVar(Form("%s_EE_totalEvents",tag.Data()),"",EE->sumEntries(),0,1e9);

  outputWs->import(*totEB);
  outputWs->import(*totEE);

  cout << EB <<endl;
  cout << EE <<endl;

  Long64_t iEntry=-1;
  
  const RooArgSet *args;

  args = EB->get(0);
  RooDataSet EB00(Form("%s_EB_0_0",tag.Data()),"",*args);
  RooDataSet EB01(Form("%s_EB_0_1",tag.Data()),"",*args);
  RooDataSet EB10(Form("%s_EB_1_0",tag.Data()),"",*args);
  RooDataSet EB11(Form("%s_EB_1_1",tag.Data()),"",*args);

  RooDataSet EE00(Form("%s_EE_0_0",tag.Data()),"",*args);
  RooDataSet EE01(Form("%s_EE_0_1",tag.Data()),"",*args);
  RooDataSet EE10(Form("%s_EE_1_0",tag.Data()),"",*args);
  RooDataSet EE11(Form("%s_EE_1_1",tag.Data()),"",*args);
  
  while( (args= EB->get(++iEntry)) ){
    if( iEntry%15000==0) cout << "Processing Entry " << iEntry <<endl;
  
    float se1 = ((RooRealVar*)args->find("sigEoE1"))->getVal();
    float pt1 = ((RooRealVar*)args->find("pt1"))->getVal();
    float etaSC1 = ((RooRealVar*)args->find("etaSC1"))->getVal();

    float se2 = ((RooRealVar*)args->find("sigEoE2"))->getVal();
    float pt2 = ((RooRealVar*)args->find("pt2"))->getVal();
    float etaSC2 = ((RooRealVar*)args->find("etaSC2"))->getVal();
    
    bool p1 = passSelection(map,se1,fabs(etaSC1),pt1);
    bool p2 = passSelection(map,se2,fabs(etaSC2),pt2);

    if( p1==0 &&  p2==0) EB00.add(*args);
    if( p1==0 &&  p2==1) EB01.add(*args);
    if( p1==1 &&  p2==0) EB10.add(*args);
    if( p1==1 &&  p2==1) EB11.add(*args);
  }

  outputWs->import(EB00);
  outputWs->import(EB01);
  outputWs->import(EB10);
  outputWs->import(EB11);


  iEntry=0;
  while( (args= EE->get(++iEntry)) ){
    if( iEntry%15000==0) cout << "Processing Entry " << iEntry <<endl;


    float se1 = ((RooRealVar*)args->find("sigEoE1"))->getVal();
    float pt1 = ((RooRealVar*)args->find("pt1"))->getVal();
    float etaSC1 = ((RooRealVar*)args->find("etaSC1"))->getVal();

    float se2 = ((RooRealVar*)args->find("sigEoE2"))->getVal();
    float pt2 = ((RooRealVar*)args->find("pt2"))->getVal();
    float etaSC2 = ((RooRealVar*)args->find("etaSC2"))->getVal();
    
    bool p1 = passSelection(map,se1,fabs(etaSC1),pt1);
    bool p2 = passSelection(map,se2,fabs(etaSC2),pt2);

    if( p1==0 &&  p2==0) EE00.add(*args);
    if( p1==0 &&  p2==1) EE01.add(*args);
    if( p1==1 &&  p2==0) EE10.add(*args);
    if( p1==1 &&  p2==1) EE11.add(*args);
  }

  outputWs->import(EE00);
  outputWs->import(EE01);
  outputWs->import(EE10);
  outputWs->import(EE11);

}

void getSB(RooWorkspace *inputWs,float offset, bool doEB,float *output){

  RooWorkspace ws("cms_hgg_spin_workspace","");

  float offsetEB = (doEB ? offset:0);
  float offsetEE = (doEB ? 0.001:offset);
  

  MakeCutWorkspaces(inputWs,&ws,"Data",offsetEB,offsetEE);
  MakeCutWorkspaces(inputWs,&ws,"Hgg125",offsetEB,offsetEE);

  if(doEB){
    MakeSignalFit(&ws,"Hgg125_EB_0_0","Hgg125_FIT_EB_0_0");
    MakeBackgroundFit(&ws,"EB_0_0","Hgg125",124.5,2,false,false);
    MakeBackgroundOnlyFit(&ws,"EB_0_0","Hgg125");
    TCanvas *c = DrawBlindFit(&ws,"EB_0_0","Hgg125",11.2);
    c->SaveAs(Form("optimization/massDist_EB_0_0_Offset_%0.5f.png",offsetEB));
    c->Delete();
  }else{
    MakeSignalFit(&ws,"Hgg125_EE_0_0","Hgg125_FIT_EE_0_0");
    MakeBackgroundFit(&ws,"EE_0_0","Hgg125",124.5,2,false,false);
    MakeBackgroundOnlyFit(&ws,"EE_0_0","Hgg125");
    TCanvas *c = DrawBlindFit(&ws,"EE_0_0","Hgg125",11.2);
    c->SaveAs(Form("optimization/massDist_EE_0_0_Offset_%0.5f.png",offsetEE));
    c->Delete();
  }
  

  double totEB  = ws.var("Hgg125_EB_totalEvents")->getVal();
  double totEE  = ws.var("Hgg125_EE_totalEvents")->getVal();

  if(doEB){
    double vEB[8];
    getSignalFitVals(&ws,"Hgg125","EB_0_0",vEB);
    
    float sigEB = vEB[0];
    float sigmaEB = vEB[6];
    
    double thisEB  = ws.data("Hgg125_EB_0_0")->sumEntries();
    
    float expSigEB = 607*thisEB/(totEB+totEE);
    
    RooAbsPdf * bkgEBf = ws.pdf("Data_Hgg125_BKGFIT_EB_0_0_bkgModel");
    RooRealVar rangeEB("rangeEB","",125-sigmaEB,125+sigmaEB);
    RooRealVar all("all","",100,180);
    float bkgEB = bkgEBf->createIntegral(rangeEB)->getVal()/bkgEBf->createIntegral(all)->getVal()*vEB[2];
    output[0] = expSigEB;
    output[1] = bkgEB;
    output[2] = thisEB;
    output[3] = vEB[2];
    output[4] = vEB[0];
    output[5] = vEB[1];
    
  }else{

    double vEE[8];
    getSignalFitVals(&ws,"Hgg125","EE_0_0",vEE);

    float sigEE = vEE[0];
    float sigmaEE = vEE[6];

    double thisEE  = ws.data("Hgg125_EE_0_0")->sumEntries();

    float expSigEE = 607*thisEE/(totEB+totEE);

    RooAbsPdf * bkgEEf = ws.pdf("Data_Hgg125_BKGFIT_EE_0_0_bkgModel");
    RooRealVar rangeEE("rangeEE","",125-sigmaEE,125+sigmaEE); 
    RooRealVar all("all","",100,180);
    float bkgEE = bkgEEf->createIntegral(rangeEE)->getVal()/bkgEEf->createIntegral(all)->getVal()*vEE[2];

    output[0] = expSigEE;
    output[1] = bkgEE;
    output[2] = thisEE;
    output[3] = vEE[2];
    output[4] = vEE[0];
    output[5] = vEE[1];
  }
}



void optimize(TString totalWSFile,bool doEB,float start,float stop,float step){

  TFile *f = new TFile(totalWSFile);
  RooWorkspace* ws = (RooWorkspace*)f->Get("cms_hgg_spin_total_workspace");

  vector<float> ssbEB,sbEB,sigEB,bkgEB,peakEB,peakErrEB;
  vector<float> ssbEE,sbEE,sigEE,bkgEE,peakEE,peakErrEE;

  float EBstart =start;
  float EBstop  =stop;
  float EBstep  =step;

  float EEstart = start;
  float EEstop  = stop;
  float EEstep  = step;

  float out[6];
  if(doEB){
    for(float EBo=EBstart;EBo<EBstop;EBo+=EBstep){
      getSB(ws,EBo,true,out);
      ssbEB.push_back(out[0]/sqrt(out[1]));
      sbEB.push_back(out[0]/out[1]);
      sigEB.push_back(out[0]);
      bkgEB.push_back(out[3]);
      peakEB.push_back(out[4]);
      peakErrEB.push_back(out[5]);
    }
  }else{
    for(float EEo=EEstart;EEo<EEstop;EEo+=EEstep){
      getSB(ws,EEo,false,out);
      ssbEE.push_back(out[0]/sqrt(out[1]));
      sbEE.push_back(out[0]/out[1]);
      sigEE.push_back(out[0]);
      bkgEE.push_back(out[3]);
      peakEE.push_back(out[4]);
      peakErrEE.push_back(out[5]);
    }
  }
  int I=0;
  cout << "offset\t\tS/sqrt(B)\tS/B\tS\tBtot\tSig\tSigErr" <<endl;
  if(doEB){
    cout << "  EB" <<endl;
    for(float EBo=EBstart;EBo<EBstop;EBo+=EBstep,I++) 
      cout << EBo << "\t\t" << ssbEB.at(I) << "\t" << sbEB.at(I) 
	   << "\t" << sigEB.at(I) << "\t" <<bkgEB.at(I) << "\t"
	   << peakEB.at(I) << "\t" << peakErrEB.at(I) <<endl;
  }else{
    cout << "\n  EE" <<endl;
    for(float EEo=EEstart;EEo<EEstop;EEo+=EEstep,I++) 
      cout << EEo << "\t\t" << ssbEE.at(I) << "\t" << sbEE.at(I) 
	   << "\t" << sigEE.at(I) << "\t" <<bkgEE.at(I) << "\t"
	   << peakEE.at(I) << "\t" << peakErrEE.at(I) <<endl;
  }

}
