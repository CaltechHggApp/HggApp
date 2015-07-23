/**

runs the whole susy hgg analysis

*/


#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <memory>
#include <cstdlib>

#include "ArgParser.hh"
#include "ReadConfig.hh"

#include "FitterNew.hpp"
#include "SMSFitterNew.hh"
#include "DataFitterNew.hh"
#include "MCBackgroundFitterNew.hh"

#include "CombinePrep.hh"

#include "ReadBinningMap.cpp"

#include "TH1D.h"
#include "TFile.h"

#include "Utils.cpp"

using namespace std;



int main(int argc,char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("ConfigFile",ArgParser::required,"path to the configuration file");
  a.addArgument("OutputFolder",ArgParser::required,"folder to save the output root files");
  
  a.addLongOption("CombineOnly",ArgParser::noArg,"run combine maker only");
  a.addLongOption("SigInjName",ArgParser::reqArg,"if CombineOnly is set, use a custom signal injection file specified as the 'data'.");
  a.addLongOption("BigSignalRegion",ArgParser::noArg,"Use a large signal region [121,129] in each category");
  a.addLongOption("nSigEffs",ArgParser::reqArg,"# of signma effectives for the signal region in each box");
  a.addLongOption("UseVariableBinning",ArgParser::noArg,"use different binning in each box depending on statistics");
  a.addLongOption("BinningMap",ArgParser::reqArg,"specify the binning map to use");


  a.addLongOption("Alternate",ArgParser::noArg,"Do the alternate analysis");

  a.addLongOption("AN239",ArgParser::noArg,"use AN13/239-like photon selection");
  a.addLongOption("highPt",ArgParser::noArg,"use 40/25 photon pT cuts");

  a.addLongOption("ForceV3Widths",ArgParser::noArg,"force the sig region widths from v3");

  a.addLongOption("HT",ArgParser::noArg,"Use HT:MET instead of razor");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  bool AN239 = a.longFlagPres("AN239");
  bool highPt = a.longFlagPres("highPt");
  bool bHT = a.longFlagPres("HT");

  bool v3Widths = a.longFlagPres("ForceV3Widths");

  bool alternate = a.longFlagPres("Alternate");

  bool externalBinning = a.longFlagPres("BinningMap");
  std::string extBinMap = "";
  if(externalBinning) extBinMap = a.getLongFlag("BinningMap");

  if(AN239 && highPt) {
    std::cout << "cannot specify both --AN239 and --highPt" << std::endl;
    return -1;
  }

  Fitter::kSelectionSet selection = Fitter::kLoose;
  if(AN239)  selection = Fitter::kAN239;
  if(highPt) selection = Fitter::kHighPt;

  string cfgFilePath = a.getArgument("ConfigFile");
  string outputFolder = a.getArgument("OutputFolder");

  bool combineOnly = a.longFlagPres("CombineOnly");

  string sigInjName="";
  if(combineOnly && a.longFlagPres("SigInjName")) sigInjName = a.getLongFlag("SigInjName");

  bool bigSigReg = a.longFlagPres("BigSignalRegion");

  float nSigEffs=2;
  if(a.longFlagPres("nSigEffs")) nSigEffs = atof( a.getLongFlag("nSigEffs").c_str() );

  cout << "Trying to make folder: " << outputFolder << "  " << mkdir(outputFolder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) << std::endl;

  ReadConfig cfg(cfgFilePath);
  cfg.printAll();
  float lumi = atof( cfg.getParameter("lumi").c_str() );

  string sm_list_string = cfg.getParameter("sm_names");
  vector<string> sm_list = cfg.tokenizeString(sm_list_string,",");

  string sms_list_string = cfg.getParameter("sms_names");
  vector<string> sms_list = cfg.tokenizeString(sms_list_string,",");

  int isMCData = atoi( cfg.getParameter("isMCData").c_str() );
  std::cout << "[INFO]: isMCData: " << isMCData << std::endl;
  string metphi = cfg.getParameter("METPhiCorrection");
  std::cout << "METPhi:  " << metphi << std::endl;

  CombinePrep combine;
  combine.setOutputFolder(outputFolder);

  std::vector<TString> catNames;
  const std::vector<TString>* cats = Fitter::getCatNames();
  catNames.resize(cats->size());
  std::copy(cats->begin(),cats->end(),catNames.begin());

  for(auto cat: catNames) std::cout << cat << std::endl;


  for(auto sm_name: sm_list) {
    string inputFileName = cfg.getParameter(sm_name+"_path");
    string strN = cfg.getParameter("N_"+sm_name);
    string strXSec = cfg.getParameter("xsec_"+sm_name);

    if(!combineOnly) {
      Fitter fitter(inputFileName,outputFolder+"/"+sm_name+".root",bHT);
      fitter.setXSec( atof(strXSec.c_str()) );
      fitter.setNTotal( atoi(strN.c_str()) );
      fitter.setLumi( lumi );
      fitter.setNSigEffs(nSigEffs);
      fitter.setSelection(selection);
      if(metphi.size()>2) fitter.setMetPhiSF( metphi.c_str() );
	
      fitter.setDoAlternateAnalysis(alternate);
      fitter.Run();
    }
    
    combine.addSMHiggsFilePath(sm_name,outputFolder+"/"+sm_name+".root");
  }

  std::map<TString,float> sigEffs; 
  //if(!combineOnly) {
  TFile SMTotFile( (outputFolder+"/SMTot.root").c_str(),"RECREATE");
    for(auto cat: catNames) {
      TH1D SMTot("SMTot_"+cat+"_mgg_hist","",3000,Fitter::minMgg,Fitter::maxMgg);
      for(auto sm_name: sm_list) {
	TFile f((outputFolder+"/"+sm_name+".root").c_str());
	SMTot.Add( (TH1D*)f.Get(cat+"_mgg_dist") );      
      }
      sigEffs[cat] = getSigEff(SMTot,125.);
      SMTot.Rebin(10);
      std::pair<float,float> width = getFWHM(SMTot);
      cout << "[INFO]: " << cat << " sigmaEff: " << sigEffs[cat] 
	   << "  wsecond-wfirst: " << width.second-width.first << std::endl;    
      SMTotFile.cd();
      SMTot.Write();
    }
    SMTotFile.Close();
    //}

    //override widths for Hbb and Zbb categories
    sigEffs["Hbb"] = 2.;
    sigEffs["Zbb"] = 2.;
    sigEffs["LowRes"] = 2.5;
    
    if(v3Widths) {
      sigEffs["HighPt"] = 1.77;
      sigEffs["HighRes"] = 1.68;
      sigEffs["LowRes"] = 2.5;
    }
    
    
  std::string dataFileName = cfg.getParameter("data_path"); 
  std::string triggerFileName = cfg.getParameter("trigger_path"); 
  std::string noiseFileName = cfg.getParameter("noise_path"); 
  if(!combineOnly) {
    DataFitter *datafitter = 0;
    if(isMCData) datafitter = new MCBackgroundFitter(dataFileName,outputFolder+"/data.root",bHT);
    else {
      std::cout << "[INFO]: ENTERING FIT..." << std::endl;
      datafitter = new DataFitter(dataFileName,outputFolder+"/data.root",bHT);
      datafitter->SetTriggerPath(triggerFileName.c_str());
      datafitter->SetNoisePath(noiseFileName.c_str());
    }
    assert(datafitter != 0);

    for(auto cat:catNames) {
      if(bigSigReg) datafitter->setSigEff(cat, 2.5);
      else datafitter->setSigEff(cat, sigEffs[cat]);
    }
    datafitter->setNSigEffs(nSigEffs);
    if(isMCData) {
      if(AN239) {
	datafitter->fixNorm("HighPt",394.6);
	datafitter->fixNorm("Hbb",3.6);
	datafitter->fixNorm("Zbb",6.8);
	datafitter->fixNorm("HighRes",934.4);
	datafitter->fixNorm("LowRes",2087.0);
      } else if(highPt) {
	datafitter->fixNorm("HighPt",654.679);
	datafitter->fixNorm("Hbb",6.68);
	datafitter->fixNorm("Zbb",6.55);
	datafitter->fixNorm("HighRes",1472.3);
	datafitter->fixNorm("LowRes",2841.8);
      } else {
	datafitter->fixNorm("HighPt",696.0);
	datafitter->fixNorm("Hbb",7.67);
	datafitter->fixNorm("Zbb",7.64);
	datafitter->fixNorm("HighRes",1553.8);
	datafitter->fixNorm("LowRes",3258.1);
      }
    }
    datafitter->setUseHT(bHT);

    datafitter->setSelection(selection);
    datafitter->setDoAlternateAnalysis(alternate);

    datafitter->Run();
    delete datafitter;
  }
  
  std::cout << "========================================" << std::endl;
  std::cout << "================Leaving Fit=============" << std::endl;
  std::cout << "========================================" << std::endl;

  if(sigInjName != "")   combine.addDataFilePath(outputFolder+"/"+sigInjName);
  else combine.addDataFilePath(outputFolder+"/data.root");

  for(auto sms_name: sms_list) {
    std::string fileName = cfg.getParameter(sms_name+"_path");
    std::string normPath = cfg.getParameter("norm_"+sms_name);

    if(!combineOnly) {
      std::cout << "========================================" << std::endl;
      std::cout << "================Entering smsFitter=============" << std::endl;
      std::cout << "========================================" << std::endl;
      SMSFitter smsFitter(fileName,outputFolder+"/"+sms_name+".root",bHT);
      smsFitter.setXSec( 1 );
      smsFitter.setLumi( lumi );
      smsFitter.setNEntriesFile( normPath );
      smsFitter.setNSigEffs(nSigEffs);

      smsFitter.setSelection(selection);
      smsFitter.setUseHT(bHT);
      smsFitter.setDoAlternateAnalysis(alternate);

      smsFitter.Run();
    }
    std::cout << "sms_name-> " << outputFolder+"/"+sms_name+".root" << std::endl;
    combine.addSMSFilePath(sms_name,outputFolder+"/"+sms_name+".root");
  }


  combine.setIsFullSIM( (TString(cfg.getParameter(sms_list[0]+"_path").c_str()).Contains("FSIM")==0) );
  combine.setCatNames(&catNames);
  combine.setSysNames(Fitter::getSysNames());

  combine.setUseVarBinning( a.longFlagPres("UseVariableBinning") );

  if(externalBinning) {
    auto binMap = ReadBinningMap(extBinMap.c_str());
    for( auto m: *binMap) {
      combine.defineExternalBinning(m.first, m.second);
    }

 }

  /*
  if(!isMCData) {
    if(bHT && AN239) {
      std::vector<int> HighPt_bins = {30,25,20,1};
      combine.defineExternalBinning("HighPt", HighPt_bins);
      std::vector<int> Hbb_bins = {8,1};
      combine.defineExternalBinning("Hbb", Hbb_bins);
      std::vector<int> Zbb_bins = {11,1};
      combine.defineExternalBinning("Zbb", Zbb_bins);
      std::vector<int> HighRes_bins = {28,21,16,1};
      combine.defineExternalBinning("HighRes", HighRes_bins);
      std::vector<int> LowRes_bins = {22,17,1};
      combine.defineExternalBinning("LowRes", LowRes_bins);      
    }else if(AN239) {
      std::vector<int> HighPt_bins = {23,16,8,3,1};
      combine.defineExternalBinning("HighPt", HighPt_bins);
      std::vector<int> Hbb_bins = {23,4,1};
      combine.defineExternalBinning("Hbb", Hbb_bins);
      std::vector<int> Zbb_bins = {23,3,1};
      combine.defineExternalBinning("Zbb", Zbb_bins);
      std::vector<int> HighRes_bins = {23,14,5,2,1};
      combine.defineExternalBinning("HighRes", HighRes_bins);
      std::vector<int> LowRes_bins = {23,14,6,3,1};
      combine.defineExternalBinning("LowRes", LowRes_bins);
    } else if(highPt) {
      std::vector<int> HighPt_bins = {23,19,8,4,2,1};
      combine.defineExternalBinning("HighPt", HighPt_bins);
      std::vector<int> Hbb_bins = {23,8,1};
      combine.defineExternalBinning("Hbb", Hbb_bins);
      std::vector<int> Zbb_bins = {23,8,1};
      combine.defineExternalBinning("Zbb", Zbb_bins);
      std::vector<int> HighRes_bins = {23,17,7,3,1};
      combine.defineExternalBinning("HighRes", HighRes_bins);
      std::vector<int> LowRes_bins = {23,17,8,4,2,1};
      combine.defineExternalBinning("LowRes", LowRes_bins);
    } else {
      std::vector<int> HighPt_bins = {23,19,8,4,2,1};
      combine.defineExternalBinning("HighPt", HighPt_bins);
      std::vector<int> Hbb_bins = {23,9,1};
      combine.defineExternalBinning("Hbb", Hbb_bins);
      std::vector<int> Zbb_bins = {23,9,1};
      combine.defineExternalBinning("Zbb", Zbb_bins);
      std::vector<int> HighRes_bins = {23,17,7,3,1};
      combine.defineExternalBinning("HighRes", HighRes_bins);
      std::vector<int> LowRes_bins = {23,17,8,4,2,1};
      combine.defineExternalBinning("LowRes", LowRes_bins);
    }
  }
  */
  combine.Make();

  system(Form("hadd -f %s/SMHiggs_SUM.root %s/ggH.root %s/vbfH.root %s/wzH.root %s/ttH.root", 
	      outputFolder.c_str(),
	      outputFolder.c_str(),
	      outputFolder.c_str(),
	      outputFolder.c_str(),
	      outputFolder.c_str()
	      ));  
  return 0;
}
