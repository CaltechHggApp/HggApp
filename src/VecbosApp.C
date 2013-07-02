//-------------------------------------------------------
// Description:
//    Routine to run Vecbos selection
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//    Maurizio Pierini
//    CERN
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>

// Vecbos includes
#include <include/Application.hh>
#include <CommonTools/include/TriggerMask.hh>
#include <include/Vecbos.hh>


// include your VECBOS analysis
#if Application == 1
#include <include/VecbosEESelection.hh>
#endif
#if Application == 2
#include <include/CandleCalib.hh>
#endif
#if Application == 3
#include <include/VecbosEESelection.hh>
#endif
#if Application == 4
#include <include/VecbosMuMuSelection.hh>
#endif
#if Application == 5
#include <include/VecbosEESelection.hh>
#endif
#if Application == 6
#include <include/VecbosMuMuSelection.hh>
#endif
#if Application == 7
#include <include/DiJet.hh>
#endif
#if Application == 8
#include <include/SUSYInclusive.hh>
#endif
#if Application == 9
#include <include/VecbosPFEESelection.hh>
#endif
#if Application == 10
#include <include/SuperClustersEESelection.hh>
#endif
//#if Application == 11
//#include <include/ElectronChecks.hh>
//#endif
#if Application == 12
#include <include/GammaPlusJet.hh>
#endif
#if Application == 13
#include <include/SUSYTau.hh>
#endif
#if Application == 14
#include <include/TopControlSample.hh>
#endif
#if Application == 15
#include <include/SUSYMultiTop.hh>
#include <include/AnalysisSelector.hh>
#endif
#if Application == 16
#include <include/CreateWJetDataset.hh>
#endif
#if Application == 17
#include <include/RazorDiPhoton.hh>
#endif
#if Application == 18
#include <include/SUSYRA.hh>
#endif
#if Application == 19
#include <include/ZMuMu.hh>
#endif
#if Application == 20
#include <include/LQ3Analysis.hh>
#endif
#if Application == 21
#include <include/VBTFLeptEff.hh>
#include <include/AnalysisSelector.hh>
#endif
#if Application == 22
#include <include/RazorLeptons.hh>
#endif
#if Application == 23
#include <include/SUSYMultiB.hh>
#include <include/AnalysisSelector.hh>
#endif
#if Application == 24
#include <include/RazorHiggsBB.hh>
#endif
#if Application == 25
#include <include/RazorBoostedTop.hh>
#endif
#if Application == 26
#include <include/SF_Filler.hh>
#include <include/AnalysisSelector.hh>
#endif
#if Application == 27
#include <include/H4b.hh>
#include <include/AnalysisSelector.hh>
#endif

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[150];
  char outFileName[150];
  char skimFileName[150];
  char json[150]="none";

  if ( argc < 3 ){
    cout << "Error at Input: please specify an input file including the list of input ROOT files" << endl; 
    cout << "Example:        ./VecbosApp list.txt output.root" << endl;
    cout << "Available options: " <<endl;
    cout << "-weight=w  weight of the MC" << endl;  
    cout << "-start=N start from event N in the chain" << endl; 
    cout << "-stop=N stop at event N in the chain" << endl; 
    cout << "-signal=N 0=W(prompt),1=Z(prompt),2=W(other),3=Z(other), 4= no mctruth" << endl; 
    cout << "--isData to run on good runs on data" << endl; 
    cout << "-json=file path of the json file" << endl;
    return 1;
  }

  // rad running options
  strcpy(inputFileName,argv[1]);
  strcpy(outFileName,argv[2]);
  strcpy(skimFileName,argv[2]);

  TChain *theChain = new TChain("ntp1");
  char Buffer[500];
  char MyRootFile[2000];  
  ifstream *inputFile = new ifstream(inputFileName);
  // get the tree with the conditions from the first file
  //  TTree *treeCond = new TTree();
  //  int nfiles=1;
  char tmpFileName[256];
  vector<string> filesToRemove;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	if(string(MyRootFile).find("eos") != std::string::npos) {
	  theChain->Add("root:/"+TString(MyRootFile));
        } else {
	  theChain->Add("rfio:"+TString(MyRootFile));
	}
        // theChain->Add("root://castorcms/"+TString(MyRootFile));
	//        theChain->Add(TString(MyRootFile));
	std::cout << "chaining " << MyRootFile << std::endl;
	//	if ( nfiles==1 ) {
	//	  TFile *firstfile = TFile::Open("root://castorcms/"+TString(MyRootFile));
	//	  treeCond = (TTree*)firstfile->Get("Conditions");
	//	}
	//        nfiles++;
      }
  }
 

  //  theChain->MakeClass("thisiswhyitcrashed");

  inputFile->close();
  delete inputFile;
  // get additional input options
  int signal = 0;
  int start = 0;
  int stop  = theChain->GetEntries();
  bool isData = false;
  float lumi = -999.;
  float xsec = -999.;
  float weight = 1.;
  for (int i=1;i<argc;i++){
    if (strncmp(argv[i],"-start",6)==0) sscanf(argv[i],"-start=%i",&start);
    if (strncmp(argv[i],"-stop",5)==0)  sscanf(argv[i],"-stop=%i",&stop);
    if (strncmp(argv[i],"-signal",7)==0)  sscanf(argv[i],"-signal=%i",&signal);
    if (strncmp(argv[i],"-weight",7)==0)  sscanf(argv[i],"-weight=%f",&weight);
    if (strncmp(argv[i],"--isData",8)==0)  isData = true;
    if (strncmp(argv[i],"-lumi",5)==0)  sscanf(argv[i],"-lumi=%f",&lumi);
    if (strncmp(argv[i],"-xsec",5)==0)  sscanf(argv[i],"-xsec=%f",&xsec);
    if (strncmp(argv[i],"-json",5)==0)  sscanf(argv[i],"-json=%s",&json);
  }

  

#if Application == 1
  std::cout << " initialising " << std::endl;
  VecbosEESelection vecbos(theChain);
  std::cout << " setup trigger selection " << std::endl;
  std::vector<string> mask;
  mask.push_back("HLT_Photon10_L1R");
  mask.push_back("HLT_Photon15_L1R");
  for(std::vector<string>::iterator it=mask.begin(); it!= mask.end();it++)
    std::cout << *it <<std::endl;
  vecbos.setRequiredTriggers(mask);
  std::cout << " setup analysis parameters" <<std::endl;
  vecbos.doBestElectronStudy(false);
  std::cout << "  skimname" <<std::endl;
  vecbos.setPrefix(skimFileName);
  std::cout << "  signal" <<std::endl;
  vecbos.setSignal(signal);
  std::cout << " starting event loop"<< std::endl;
  vecbos.Loop();
  std::cout << " preparing plots" <<std::endl;
  vecbos.displayEfficiencies();
#endif

#if Application == 2
  CandleCalib vecbos(theChain);
  //  TriggerMask mask(treeCond);
  //  mask.requireTrigger("HLT_Mu9");
  //  mask.requireTrigger("HLT_Mu11");
  //  mask.requireTrigger("HLT_Mu15");
  //  mask.requireTrigger("HLT_DoubleMu3");
  //  std::vector<int> requiredTriggers = mask.getBits();
  //  vecbos.requireTrigger(requiredTriggers);
  vecbos.Loop(string(outFileName), start, stop);  
#endif

#if Application == 3

  // for electron analysis:
  VecbosEESelection vecbos(theChain);
  std::vector<std::string> mask;
  if(!isData) { // for Winter10 MC only
    mask.push_back("HLT_Ele10_SW_L1R_v2");
    mask.push_back("HLT_Ele17_SW_L1R_v2");
    mask.push_back("HLT_Ele17_SW_L1R_v2");
    mask.push_back("HLT_Ele17_SW_Isol_L1R_v2");
    mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    mask.push_back("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2");
  } else { // 2010 Run data
    std::cout << "Applying triggers for data" << std::endl;
    mask.push_back("HLT_Photon10_L1R");
    mask.push_back("HLT_Photon15_L1R");
    mask.push_back("HLT_Photon15_Cleaned_L1R");
    mask.push_back("HLT_Ele10_LW_L1R");
    mask.push_back("HLT_Ele15_SW_L1R");
    mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
    mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
    mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
    mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
    mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
  }
  vecbos.setRequiredTriggers(mask);
  vecbos.doBestElectronStudy(false);
  vecbos.doBTagEfficiencyStudy(false);
  vecbos.setPrefix(skimFileName);
  vecbos.setSignal(signal);
  vecbos.Loop();
  vecbos.displayEfficiencies();

  return 0;

#endif

#if Application == 4  
  // for muon analysis:
  VecbosMuMuSelection vecbos(theChain);
  std::vector<std::string> mask;
  mask.push_back("HLT_Mu9");
  mask.push_back("HLT_Mu15_v1");
  vecbos.setRequiredTriggers(mask);
  vecbos.setPrefix(skimFileName);
  vecbos.setSignal(signal);
  vecbos.Loop();
  vecbos.displayEfficiencies();
  return 0;
#endif

#if Application == 5
#endif

#if Application == 6
  
#endif

#if Application == 7
  if(isData) {
    DiJet vecbos(theChain, true, true);
    vecbos.Loop(string(outFileName), start, stop);  
  } else {
    DiJet vecbos(theChain, false, false);
    vecbos.Loop(string(outFileName), start, stop);  
  }
#endif

#if Application == 8  
  if(isData) {
    vecbos.Loop(string(outFileName), start, stop);  
  } else {
    Razor vecbos(theChain, string(json), false, false);
    vecbos.Loop(string(outFileName), start, stop);  
  }
#endif

#if Application == 9

  VecbosPFEESelection vecbos(theChain);
  std::vector<std::string> mask;
  mask.push_back("HLT_Photon10_L1R");
  vecbos.setRequiredTriggers(mask);
  vecbos.setPrefix(skimFileName);
  vecbos.setSignal(signal);
  vecbos.Loop();
  vecbos.displayEfficiencies();

#endif

#if Application == 10

  // for electron analysis:
  SuperClustersEESelection vecbos(theChain);
  std::vector<std::string> mask;
  mask.push_back("HLT_Photon10_L1R");
  vecbos.setRequiredTriggers(mask);
  vecbos.requireTrigger(requiredTriggers);
  vecbos.setPrefix(skimFileName);
  vecbos.setSignal(signal);
  vecbos.Loop();
  vecbos.displayEfficiencies();
#endif
 //  **********NATASHA MODIFIED 6/21/2010*****************

#if Application == 12
  if(isData) {
    GammaPlusJet vecbos(theChain, true, true);
    vecbos.SetConditions(treeCond);
    vecbos.Loop(string(outFileName), start, stop);  
  } else {
    GammaPlusJet vecbos(theChain, false, false);
    if(lumi > 0. && xsec > 0.) {
      vecbos.SetLuminosity(lumi);
      vecbos.SetXsection(xsec);
    }
    vecbos.SetConditions(treeCond);
    vecbos.Loop(string(outFileName), start, stop);  
  }
#endif

#if Application == 13  
  std::vector<std::string> mask;
  // if(!isData) {
  //   mask.push_back("HLT_HT200");
  // } else {
  if (isData) {
    mask.push_back("HLT_R014_MR150");
    mask.push_back("HLT_R020_MR150");
    mask.push_back("HLT_R025_MR150");
    mask.push_back("HLT_R020_MR500");
    mask.push_back("HLT_R020_MR550");
    mask.push_back("HLT_R025_MR400");
    mask.push_back("HLT_R025_MR450");
    mask.push_back("HLT_R033_MR300");
    mask.push_back("HLT_R033_MR350");
    mask.push_back("HLT_R038_MR200");
    mask.push_back("HLT_R038_MR250");
    mask.push_back("HLT_Mu8_R005_MR200");
    mask.push_back("HLT_Mu8_R020_MR200");
    mask.push_back("HLT_Mu8_R025_MR200");
    mask.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200");
    mask.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200");
    mask.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200");
  }
  if(isData) {
    SUSYInclusive vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);  
  } else {
    SUSYTau vecbos(theChain, false, false);
    // vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);  
  }
#endif

#if Application == 14

  // top control sample
  TopControlSample vecbos(theChain);
  std::vector<std::string> mask;
  mask.push_back("HLT_Photon10_L1R");
  mask.push_back("HLT_Photon15_L1R");
  mask.push_back("HLT_Photon15_Cleaned_L1R");
  mask.push_back("HLT_Ele10_LW_L1R");
  mask.push_back("HLT_Ele15_SW_L1R");
  mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
  mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
  mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
  mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
  mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");

  mask.push_back("HLT_Mu9");
  mask.push_back("HLT_Mu15_v1");

  // MC triggers
//   mask.push_back("HLT_Ele15_LW_L1R");
//   mask.push_back("HLT_Mu9");
  vecbos.setRequiredTriggers(mask);
  vecbos.setPrefix(skimFileName);
  vecbos.setSignal(signal);
  vecbos.Loop();
  vecbos.displayEfficiencies();

#endif


#if Application == 15
  std::vector<std::string> mask;
  if(isData) {
    if(AnalysisSelector == 1 || AnalysisSelector == 2){

      //      mask.push_back("HLT_Mu9");
      //      mask.push_back("HLT_Mu11");
      //      mask.push_back("HLT_Mu15");                                                                        

      //Single Mu
      // mask.push_back("1-163261:HLT_Mu15_v2");
//       mask.push_back("163262-164237:HLT_Mu24_v");
//       mask.push_back("165085-999999:HLT_Mu30_v");
      //Iso Mu
           mask.push_back("1-163261:HLT_Mu15_v2");
           mask.push_back("163262-167043:HLT_IsoMu17_v");
           mask.push_back("167044-167913:HLT_IsoMu17_eta2p1_v");
           mask.push_back("170053-172949:HLT_IsoMu20_v");
      
    }else if(AnalysisSelector == 3 || AnalysisSelector == 4){
      mask.push_back("HLT_Ele10_LW_L1R");                                                                                   
      mask.push_back("HLT_Ele15_SW_L1R");
      mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");                                                               
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    }
    SUSYMultiTop vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);  
  } else {
    SUSYMultiTop vecbos(theChain, false, false);
    vecbos.SetWeight(double(weight));
    vecbos.Loop(string(outFileName), start, stop);  
  }
#endif

#if Application == 23
  std::vector<std::string> mask;
  if(isData) {
    if(AnalysisSelector == 1 || AnalysisSelector == 2){
      mask.push_back("1-163261:HLT_Mu15_v2");                                                                                                      
      mask.push_back("163262-164237:HLT_Mu24_v"); 
      mask.push_back("165085-999999:HLT_Mu30_v"); 
      
    }else if(AnalysisSelector == 3 || AnalysisSelector == 4){
      mask.push_back("HLT_Ele10_LW_L1R");
      mask.push_back("HLT_Ele15_SW_L1R");
      mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    }
    SUSYMultiB vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    SUSYMultiB vecbos(theChain, false, false);
    vecbos.SetWeight(double(weight));
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif


#if Application == 16
  if(isData == true)
    {
      CreateWJetDataset vecbos(theChain, true, true, -1);
      std::vector<std::string> mask;
      mask.push_back("HLT_Mu15");
      mask.push_back("HLT_Mu15_v1");
      vecbos.setRequiredTriggers(mask);
      vecbos.Loop(string(outFileName), start, stop); 
    }
  else
    {
      int SelectZMuMu = -1;
      if(signal == 1)   // Z(prompt)
	SelectZMuMu = 1;
      if(signal == 3)   // Z(other)
	SelectZMuMu = 0;
      
      CreateWJetDataset vecbos(theChain, false, false, SelectZMuMu);
      std::vector<std::string> mask;
      mask.push_back("HLT_Mu15");
      mask.push_back("HLT_Mu15_v1");
      vecbos.setRequiredTriggers(mask);
      vecbos.Loop(string(outFileName), start, stop);
    }
#endif

#if Application == 17
  if(isData) {
    RazorDiPhoton vecbos(theChain, string(json), true, true);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    RazorDiPhoton vecbos(theChain,  string(json), false, false);
    vecbos.SetWeight(double(weight));
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif

#if Application == 18
  std::vector<std::string> mask;
  if(isData) {
    if(RA4Selector == 1){
      mask.push_back("HLT_Mu11");
      mask.push_back("HLT_Mu5_HT70U_v3");
    }else if(RA4Selector == 2){
      mask.push_back("HLT_Ele15_SW_L1R");
      mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
      mask.push_back("HLT_Ele10_SW_HT70U_L1R_v1");
      mask.push_back("HLT_Ele10_SW_HT70U_L1R_v2");
    }
    SUSYRA vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    SUSYRA vecbos(theChain, false, false);
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif

#if Application == 19
  std::vector<std::string> mask;
  if(isData) {
    mask.push_back("HLT_Mu9");
    mask.push_back("HLT_Mu11");
    mask.push_back("HLT_Mu15");
    mask.push_back("HLT_Mu17");
    ZMuMu vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    ZMuMu vecbos(theChain, false, false);
    vecbos.Loop(string(outFileName), start, stop);
  }

#endif

#if Application == 20   // LQ3 Analysis
  std::vector<std::string> mask;
  if(isData == true)
  {
     LQ3Analysis Analysis(theChain, true, true);
     Analysis.setRequiredTriggers(mask);
     Analysis.Loop(string(outFileName), start, stop);
  }
  else
  {
     LQ3Analysis Analysis(theChain, false, true);
     Analysis.Loop(string(outFileName), start, stop);
  }
#endif

#if Application == 21   // VBTF Lepton Efficiencies                                                                                                                     
  std::vector<std::string> mask;
  if(isData) {
    if(AnalysisSelector == 1 || AnalysisSelector == 2){
      mask.push_back("HLT_Mu9");
      mask.push_back("HLT_Mu11");
      mask.push_back("HLT_Mu15");
    }else if(AnalysisSelector == 3 || AnalysisSelector == 4){
      mask.push_back("HLT_Ele10_LW_L1R");
      mask.push_back("HLT_Ele15_SW_L1R");
      mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    }
    VBTFLeptEff vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    VBTFLeptEff vecbos(theChain, false, false);
    vecbos.Loop(string(outFileName), start, stop);
  }

#endif

#if Application == 22
  if(isData) {
    RazorLeptons vecbos(theChain, string(json), true, true);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    RazorLeptons vecbos(theChain, string(json), false, false);
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif

#if Application == 24
  if(isData) {
    RazorHiggsBB vecbos(theChain, string(json), true, true);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    RazorHiggsBB vecbos(theChain, string(json), false, false);
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif

#if Application == 25
  if(isData) {
    RazorBoostedTop vecbos(theChain, string(json), true, true);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    int model=-1;
    const int nModels=7;
    const char *models[nModels] = {"mSUGRA","T1bbbb","T2bb","T2tt","T1tttt","T2","T1"};
    int modNameLen[nModels] = {6,6,4,4,6,2,2};
    
    for(int i=0;i<nModels;i++){
      if(strstr(inputFileName,models[i])){
	model=i; break;
      }
    }
    RazorBoostedTop vecbos(theChain, string(json), false, false,model);
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif

#if Application == 26
  std::vector<std::string> mask;
  if(isData) {
    if(AnalysisSelector == 1 || AnalysisSelector == 2){
      mask.push_back("1-163261:HLT_Mu15_v2");
      mask.push_back("163262-164237:HLT_Mu24_v");
      mask.push_back("165085-999999:HLT_Mu30_v");

    }else if(AnalysisSelector == 3 || AnalysisSelector == 4){
      mask.push_back("HLT_Ele10_LW_L1R");
      mask.push_back("HLT_Ele15_SW_L1R");
      mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    }
    SF_Filler vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    SF_Filler vecbos(theChain, false, false);
    vecbos.SetWeight(double(weight));
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif

#if Application == 27
  std::vector<std::string> mask;
  if(isData) {
    if(AnalysisSelector == 1 || AnalysisSelector == 2){
      mask.push_back("1-163261:HLT_Mu15_v2");
      mask.push_back("163262-164237:HLT_Mu24_v");
      mask.push_back("165085-999999:HLT_Mu30_v");

    }else if(AnalysisSelector == 3 || AnalysisSelector == 4){
      mask.push_back("HLT_Ele10_LW_L1R");
      mask.push_back("HLT_Ele15_SW_L1R");
      mask.push_back("HLT_Ele15_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_CaloEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TightEleId_L1R");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2");
      mask.push_back("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    }
    H4b vecbos(theChain, true, true);
    vecbos.setRequiredTriggers(mask);
    vecbos.Loop(string(outFileName), start, stop);
  } else {
    H4b vecbos(theChain, false, false);
    vecbos.SetWeight(double(weight));
    vecbos.Loop(string(outFileName), start, stop);
  }
#endif


  system("rm thisiswhyitcrashed*");

  return 0;

}
