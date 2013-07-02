#include "ArgParser.hh"
#include "MakeSpinPlots.h"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  
  ArgParser a(argc,argv);
  a.addArgument("InputWorkspace",ArgParser::required,"input workspace");
  a.addArgument("Lumi",ArgParser::required,"Luminosity");
  a.addArgument("OutputPath",ArgParser::required,"output directory");
  a.addArgument("OutputTag",ArgParser::required,"tag for the output plots");
  a.addLongOption("PrintOnly",ArgParser::noArg,"Only print the yields, don't make plots");
  a.addLongOption("AllMC",ArgParser::noArg,"make plots for all MC samples in the workspace (default: just Hgg125)");
  a.addLongOption("workspace",ArgParser::required,"name of the workspace to process (default: cms_hgg_spin_workspace)");
  a.addLongOption("SMName",ArgParser::required,"name of the sample to draw as the SM Higgs (default: none)");

  string ret;
  if(a.process(ret) != 0){
    a.printOptions(argv[0]);
    return 0;
  }

  string inputWS = a.getArgument("InputWorkspace");
  float lumi = atof(a.getArgument("Lumi").c_str());
  string bp  = a.getArgument("OutputPath");
  string tag = a.getArgument("OutputTag");
  bool pOnly = a.longFlagPres("PrintOnly");
  bool all = a.longFlagPres("AllMC");
  MakeSpinPlots *msp;
  if(a.longFlagPres("workspace")){
    msp = new MakeSpinPlots(inputWS,tag,a.getLongFlag("workspace"));
  }else{
    msp = new MakeSpinPlots(inputWS,tag);    
  }
  msp->setLumi(lumi);
  msp->setBasePath(bp);

  if(a.longFlagPres("SMName")) msp->setSMName(a.getLongFlag("SMName"));
  if(!pOnly){
    if(all) msp->runAll();
    else msp->runAll("Hgg125");
  }

  msp->printAll();
  delete msp;
}
