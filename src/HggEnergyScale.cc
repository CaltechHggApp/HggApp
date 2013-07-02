#include "HggEnergyScale.hh"
#include <iostream>
using namespace std;

#define debugEnergyScale 0

HggEnergyScale::HggEnergyScale(std::string path){
  valid= false; // only s

  ReadConfig reader;
  try{
    reader.read(path);
  }catch(std::exception &e){
    std::cout << "HggEnergyScale" <<std::endl;
    std::cout << path <<std::endl;
    throw  e;
  }


  configNames.push_back("EBlowEtaBadDeltaE");
  configNames.push_back("EBlowEtaGoldDeltaE");
  configNames.push_back("EBhiEtaBadDeltaE");
  configNames.push_back("EBhiEtaGoldDeltaE");
  configNames.push_back("EElowEtaBadDeltaE");
  configNames.push_back("EElowEtaGoldDeltaE");
  configNames.push_back("EEhiEtaBadDeltaE");
  configNames.push_back("EEhiEtaGoldDeltaE");

  configNames.push_back("EBlowEtaBadDeltaE_Err");
  configNames.push_back("EBlowEtaGoldDeltaE_Err");
  configNames.push_back("EBhiEtaBadDeltaE_Err");
  configNames.push_back("EBhiEtaGoldDeltaE_Err");
  configNames.push_back("EElowEtaBadDeltaE_Err");
  configNames.push_back("EElowEtaGoldDeltaE_Err");
  configNames.push_back("EEhiEtaBadDeltaE_Err");
  configNames.push_back("EEhiEtaGoldDeltaE_Err");

  configNames.push_back("EBlowEtaBadDeltaE_ScaleErr");
  configNames.push_back("EBlowEtaGoldDeltaE_ScaleErr");
  configNames.push_back("EBhiEtaBadDeltaE_ScaleErr");
  configNames.push_back("EBhiEtaGoldDeltaE_ScaleErr");
  configNames.push_back("EElowEtaBadDeltaE_ScaleErr");
  configNames.push_back("EElowEtaGoldDeltaE_ScaleErr");
  configNames.push_back("EEhiEtaBadDeltaE_ScaleErr");
  configNames.push_back("EEhiEtaGoldDeltaE_ScaleErr");

  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);

  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);

  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);
  highR9.push_back(false);
  highR9.push_back(true);

  minEta.push_back(0.);
  minEta.push_back(0.);
  minEta.push_back(1.);
  minEta.push_back(1.);
  minEta.push_back(1.48);
  minEta.push_back(1.48);
  minEta.push_back(2.);
  minEta.push_back(2.);

  minEta.push_back(0.);
  minEta.push_back(0.);
  minEta.push_back(1.);
  minEta.push_back(1.);
  minEta.push_back(1.48);
  minEta.push_back(1.48);
  minEta.push_back(2.);
  minEta.push_back(2.);

  minEta.push_back(0.);
  minEta.push_back(0.);
  minEta.push_back(1.);
  minEta.push_back(1.);
  minEta.push_back(1.48);
  minEta.push_back(1.48);
  minEta.push_back(2.);
  minEta.push_back(2.);

  maxEta.push_back(1.);
  maxEta.push_back(1.);
  maxEta.push_back(1.48);
  maxEta.push_back(1.48);
  maxEta.push_back(2.);
  maxEta.push_back(2.);
  maxEta.push_back(3.);
  maxEta.push_back(3.);

  maxEta.push_back(1.);
  maxEta.push_back(1.);
  maxEta.push_back(1.48);
  maxEta.push_back(1.48);
  maxEta.push_back(2.);
  maxEta.push_back(2.);
  maxEta.push_back(3.);
  maxEta.push_back(3.);

  maxEta.push_back(1.);
  maxEta.push_back(1.);
  maxEta.push_back(1.48);
  maxEta.push_back(1.48);
  maxEta.push_back(2.);
  maxEta.push_back(2.);
  maxEta.push_back(3.);
  maxEta.push_back(3.);

  nRegions = configNames.size();

  if(debugEnergyScale) cout << "reader status: " << reader.is_init() << endl;
  if(!reader.is_init()) return;

  //std::vector<TString> strings;
  //for(int i=0;i<nRegions;i++) strings.push_back( (TString)reader.getParameter(configNames[i]) );

  char runString[400];
  strcpy(runString,reader.getParameter("RunUpper").c_str()); //need a non-const char*
  //tokenize the string
  if(debugEnergyScale) cout << runString << endl;
  char *rs = strtok(runString,",");
  while(rs){
    if(debugEnergyScale) cout << rs << endl;
    runs.push_back(atoi(rs));
    rs = strtok(NULL,",");
  }


  char valString[400];
  for(int i=0;i<nRegions;i++){
    strcpy(valString,reader.getParameter(configNames[i]).c_str());
    energyScales[configNames[i]] = std::vector<float>();
    char *vs = strtok(valString,",");
    while(vs){
      if(debugEnergyScale) cout << "Region " << i << ": " << vs << endl;
      energyScales[configNames[i]].push_back(atof(vs));
      vs = strtok(NULL,",");
    }
  }

  valid = true;
}

int HggEnergyScale::getRunIndex(int run){
  if(debugEnergyScale) cout << "getRunIndex" << endl;
  int runIndex;
  for(runIndex=0;runIndex<runs.size();runIndex++){
    if(run <= runs.at(runIndex)) break;
  }
  if(runIndex == runs.size()) runIndex--;
  if(debugEnergyScale) cout << runIndex << endl;
  if(debugEnergyScale) cout << "run index: " << runIndex << endl;
  return runIndex;
}

std::pair<float,float> HggEnergyScale::getDEoE(VecbosPho pho, int run){
  if(debugEnergyScale) cout << "getDEoE" << endl;
  if(!valid) return std::pair<float,float>(0,0);
  int runIndex = getRunIndex(run);
  int selectRegion = getCategory(pho);

  if(selectRegion == -1) return std::pair<float,float>(0,0);

  string regionName = configNames[selectRegion];
  string regionErrName = regionName;  regionErrName.append("_Err");
  
  if(runIndex<0 || runIndex >= energyScales[regionName].size() ) return std::pair<float,float>(1,0);

  if(debugEnergyScale) cout << energyScales[regionName].at(runIndex) << endl;
  return std::pair<float,float>(energyScales[regionName].at(runIndex),
				energyScales[regionErrName].at(runIndex));
}

float HggEnergyScale::getMCScaleErr(VecbosPho pho, int run){
  if(debugEnergyScale) cout << "getMCScaleErr" << endl;
  int runIndex = getRunIndex(run);
  int selectRegion = getCategory(pho);

  if(selectRegion == -1) return -999;

  string regionName = configNames[selectRegion];  regionName.append("_ScaleErr");
  return energyScales[regionName].at(runIndex);
}

float HggEnergyScale::getCategory(VecbosPho pho){
  if(debugEnergyScale) cout << "getCategory" << endl;
  int selectRegion=-1;
  //cout << pho.eta << "  " << pho.SC.r9 << endl;
  for(int iReg = 0; iReg< nRegions; iReg++){
    //cout << ">> " << minEta[iReg] << "  " << maxEta[iReg] << "  " << r9Cut << "  " << highR9[iReg] << endl; 
    if( fabs(pho.SC.eta) >= minEta[iReg] 
	&& fabs(pho.SC.eta) <maxEta[iReg] 
	&& ((pho.SC.r9 > r9Cut) == highR9[iReg]) ){
      selectRegion = iReg;
      break;
    }
  }
  if(debugEnergyScale) cout << selectRegion << endl;
  return selectRegion;
}

