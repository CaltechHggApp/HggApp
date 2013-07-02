#include "ArgParser.hh"
#include "MakeSpinWorkspace.h"
#include "ReadConfig.hh"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  ArgParser a(argc,argv);
  a.addArgument("WorkspaceFile",ArgParser::required,"path to the workspace file");
  a.addArgument("ConfigFile",  ArgParser::required,"path to the config file giving the data and MC options");
  a.addArgument("ConfigOption",ArgParser::required,"which heading in the config file");
  a.addLongOption("SelectionMap",ArgParser::reqArg,"Which selection map to use");
  a.addLongOption("DataEfficiencyMap",ArgParser::reqArg,"Which map to use for efficiency correction for data (default: no efficiency correction)");
  a.addLongOption("MCEfficiencyMap",ArgParser::reqArg,"Which map to use for efficiency correction for MC (default: no efficiency correction)");
  a.addLongOption("noCiC",ArgParser::noArg,"specify to disable CiC selection (default: on)");
  a.addLongOption("useR9",ArgParser::noArg,"use r9 categories (default: off)");
  a.addLongOption("tightPt",ArgParser::noArg,"use tight pt/m cuts (default: off)");
  a.addLongOption("mixMC",ArgParser::required,"Specify two MC samples to mix, comma-separated (e.g.: --mixMC=Hgg125,RSG125)");
  a.addLongOption("fractions",ArgParser::required,"Specify the fractions of the MC samples specified in --mixMC, comma-separated (e.g.: --fractions=0.1,0.2,0.8)");
  a.addLongOption("useUncorrectedMass",ArgParser::noArg,"Use the uncorrected (no regression, scale or smear) mass");
  a.addLongOption("setMassRange",ArgParser::reqArg,"set the mass range to include (default 100-180)");
  a.addLongOption("setLeadPtCut",ArgParser::reqArg,"set the pT cut of the leading photon (default: 32)");
  a.addLongOption("setTrailingPtCut",ArgParser::reqArg,"set the pT cut of the trailing photon (default: 24)");
  a.addLongOption("useHelicityFrame",ArgParser::noArg,"use the helicity frame to measure cos(theta) (default: collins-sopper frame)");
  a.addLongOption("useAsymmCosTheta",ArgParser::noArg,"use the asymmetric (non-absolute value) cos(theta)");
  a.addLongOption("noData",ArgParser::noArg,"don't run data (useful for adding additional MC samples)");
  a.addLongOption("catFromTree",ArgParser::noArg,"take the event category from input tree (globe only, default: no)");
  a.addLongOption("optimization",ArgParser::noArg,"workspace is for optimization, saves more variables");
  a.addLongOption("TwoEBCats",ArgParser::noArg,"split the EB into two categories at eta=0.8");
  a.addLongOption("VetoInnerEE",ArgParser::noArg,"Veto the inner region of the EE");
  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  string wsFile = a.getArgument("WorkspaceFile");

  string cfgFile = a.getArgument("ConfigFile");
  string cfgOption = a.getArgument("ConfigOption");


  //
  // READ CONFIGURATION FILE
  //
  ReadConfig cfgReader(cfgFile,ReadConfig::kSection);
  string data = cfgReader.getParameter("data",cfgOption);
  string mcList = cfgReader.getParameter("mcList",cfgOption);
  cout << mcList <<endl;
  vector<string> mcListVec= cfgReader.tokenizeString(mcList,",");
  cout << mcListVec.size() <<endl;
  int runMin = atoi(cfgReader.getParameter("runMin",cfgOption).c_str());
  int runMax = atoi(cfgReader.getParameter("runMax",cfgOption).c_str());

  string KFactorFile = cfgReader.getParameter("KFactorFile",cfgOption);
  string RescaleFile = cfgReader.getParameter("RescaleFile",cfgOption);

  bool isGlobe = (cfgReader.getParameter("isGlobe",cfgOption).compare("yes")==0);
  bool isList = (cfgReader.getParameter("isList",cfgOption).compare("yes")==0);

  float lumi = atof(cfgReader.getParameter("lumi",cfgOption).c_str());

  if(lumi<=0){
    std::cout << "LUMI NOT SPECIFIED FOR " << cfgOption << std::endl
	      << "ABORTING" <<std::endl;
    return 1;
  }

  if(isGlobe){
    std::cout << "\n\n RUNNING ON GLOBE NTUPLES\n\n" << std::endl;
  }
  //
  // ------
  //


  int selectionMap=7;
  if(a.longFlagPres("SelectionMap")) selectionMap = atoi(a.getLongFlag("SelectionMap").c_str());
  bool requireCiC=true;
  if(a.longFlagPres("noCiC")) requireCiC=false;
  bool useR9 = a.longFlagPres("useR9");
  bool tightPt = a.longFlagPres("tightPt");
  TString effMap_data="";
  TString effMap_mc="";
  if(a.longFlagPres("DataEfficiencyMap")) effMap_data = a.getLongFlag("DataEfficiencyMap").c_str();
  if(a.longFlagPres("MCEfficiencyMap")) effMap_mc = a.getLongFlag("MCEfficiencyMap").c_str();


  MakeSpinWorkspace msw(wsFile);

  if(a.longFlagPres("mixMC") && a.longFlagPres("fractions")){
    msw.setMixDatasets();
    std::vector<string> mc = ReadConfig::tokenizeString(a.getLongFlag("mixMC"),",");
    std::vector<string> fst = ReadConfig::tokenizeString(a.getLongFlag("fractions"),",");
    cout << mc.at(0) << endl << mc.at(1) <<endl;
    if(mc.size() != 2){
      cout << "\n\n invalid argument to flag --mixMC\n" <<std::endl;
      a.printOptions(argv[0]);
      return -1;
    }
    if(fst.size()==0){
      cout << "\n\n invalid argument to flag --fractions\n" <<std::endl;
      a.printOptions(argv[0]);
      return -1;
    }

    for(int i=0;i<fst.size();i++){
      msw.getMixer()->scheduleMix(mc.at(0).c_str(),mc.at(1).c_str(),atof(fst.at(i).c_str()));
    }
  }
  if(a.longFlagPres("useUncorrectedMass")) msw.setUseUncorrMass();

  cout << "Data:    " << data <<endl;


  if(!a.longFlagPres("noData"))   msw.addFile(data,"Data",true,-1,isList);
  for(vector<string>::const_iterator mcIt = mcListVec.begin();
      mcIt != mcListVec.end();
      mcIt++){
    string mcName = *mcIt;
    string filePath = cfgReader.getParameter(mcName,cfgOption);
    string NgenSt   = cfgReader.getParameter(Form("N_%s",mcName.c_str()),cfgOption);
    int Ngen = -1;
    if(NgenSt.compare("")!=0) Ngen = atoi(NgenSt.c_str());
    cout << mcName << ":    " << filePath <<endl;
    msw.addFile(filePath,mcName,false,Ngen,isList);
  }
  
  msw.setIsGlobe(isGlobe);
  if(a.longFlagPres("catFromTree")) msw.setTakeCatFromTree();
  if(a.longFlagPres("optimization")) msw.setOptimization();
  if(a.longFlagPres("TwoEBCats")) msw.setTwoEBCats();
  if(a.longFlagPres("VetoInnerEE")) msw.setVetoInnerEE();


  msw.setLumi(lumi);

  msw.setRequireCiC(requireCiC);
  msw.setSelectionMap(selectionMap);
  msw.setRunRange(runMin,runMax);
  std::cout << "Data Efficiency Correction: " << effMap_data << std::endl
	    << "MC   Efficiency Correction: " << effMap_mc   << std::endl;
  msw.setEfficiencyCorrectionFile(effMap_data,effMap_mc);
  msw.setUseR9(useR9);
  msw.setTightPt(tightPt);
  msw.setKFactorFile(KFactorFile.c_str());
  msw.setRescaleFile(RescaleFile.c_str());

  if(a.longFlagPres("setMassRange")){
    std::vector<string> range = ReadConfig::tokenizeString(a.getLongFlag("setMassRange"),"-");
    if(range.size() != 2){
      cout << "\n\n Invalid argument to flag --setMassRange" <<endl;
      a.printOptions(argv[0]);
      return -1;      
    }
    msw.setMassRange( atof(range.at(0).c_str()),atof(range.at(1).c_str()) );
  }

  msw.setUseHelicityFrame(a.longFlagPres("useHelicityFrame"));
  msw.setUseAbsCosTheta(!a.longFlagPres("useAsymmCosTheta"));

  msw.MakeWorkspace();

  return 0;

}
