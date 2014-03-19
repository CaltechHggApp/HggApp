/**

runs the whole susy hgg analysis

*/


#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "ArgParser.hh"
#include "ReadConfig.hh"

#include "FitterNew.hpp"
#include "SMSFitterNew.hh"
#include "DataFitterNew.hh"

#include "CombinePrep.hh"


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
  a.addLongOption("UseVariableBinning",ArgParser::noArg,"use different binning in each box depending on statistics");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  string cfgFilePath = a.getArgument("ConfigFile");
  string outputFolder = a.getArgument("OutputFolder");

  bool combineOnly = a.longFlagPres("CombineOnly");

  string sigInjName="";
  if(combineOnly && a.longFlagPres("SigInjName")) sigInjName = a.getLongFlag("SigInjName");

  bool bigSigReg = a.longFlagPres("BigSignalRegion");

  cout << "Trying to make folder: " << outputFolder << "  " << mkdir(outputFolder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) << std::endl;

  ReadConfig cfg(cfgFilePath);
  cfg.printAll();
  float lumi = atof( cfg.getParameter("lumi").c_str() );

  string sm_list_string = cfg.getParameter("sm_names");
  vector<string> sm_list = cfg.tokenizeString(sm_list_string,",");

  string sms_list_string = cfg.getParameter("sms_names");
  vector<string> sms_list = cfg.tokenizeString(sms_list_string,",");

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
      Fitter fitter(inputFileName,outputFolder+"/"+sm_name+".root");
      fitter.setXSec( atof(strXSec.c_str()) );
      fitter.setNTotal( atoi(strN.c_str()) );
      fitter.setLumi( lumi );
      
      fitter.Run();
    }

    combine.addSMHiggsFilePath(sm_name,outputFolder+"/"+sm_name+".root");
  }

  std::map<TString,float> sigEffs; 
  //if(!combineOnly) {
  TFile SMTotFile( (outputFolder+"/SMTot.root").c_str(),"RECREATE");
    for(auto cat: catNames) {
      TH1D SMTot("SMTot_"+cat+"_mgg_hist","",3000,110,140);
      for(auto sm_name: sm_list) {
	TFile f((outputFolder+"/"+sm_name+".root").c_str());
	SMTot.Add( (TH1D*)f.Get(cat+"_mgg_dist") );      
      }
      sigEffs[cat] = getSigEff(SMTot,125.);
      SMTot.Rebin(10);
      std::pair<float,float> width = getFWHM(SMTot);
      cout << cat << "  " << sigEffs[cat] << "  " << width.second-width.first << std::endl;    
      SMTotFile.cd();
      SMTot.Write();
    }
    SMTotFile.Close();
    //}

  std::string dataFileName = cfg.getParameter("data_path"); 
  if(!combineOnly) {
    DataFitter datafitter(dataFileName,outputFolder+"/data.root");
    for(auto cat:catNames) {
      if(bigSigReg) datafitter.setSigEff(cat, 2.0);
      else datafitter.setSigEff(cat, sigEffs[cat]);
    }
    datafitter.Run();

  }
  if(sigInjName != "")   combine.addDataFilePath(outputFolder+"/"+sigInjName);
  else combine.addDataFilePath(outputFolder+"/data.root");

  for(auto sms_name: sms_list) {
    std::string fileName = cfg.getParameter(sms_name+"_path");
    std::string normPath = cfg.getParameter("norm_"+sms_name);

    if(!combineOnly) {
      SMSFitter smsFitter(fileName,outputFolder+"/"+sms_name+".root");
      smsFitter.setXSec( 1 );
      smsFitter.setLumi( lumi );
      smsFitter.setNEntriesFile( normPath );
      
      smsFitter.Run();
    }
    combine.addSMSFilePath(sms_name,outputFolder+"/"+sms_name+".root");
  }


  combine.setCatNames(&catNames);
  combine.setSysNames(Fitter::getSysNames());

  combine.setUseVarBinning( a.longFlagPres("UseVariableBinning") );

  combine.Make();
  
  return 0;
}
