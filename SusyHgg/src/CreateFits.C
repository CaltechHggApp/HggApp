
#include <iostream>
#include <string>
#include "exception"
#include <stdexcept>


#include "RooWorkspace.h"

#include "TFile.h"

#include "../include/ArgParser.hh"

#include "RunFits.C"
#include "Fitter.hpp"

int main(int argc,char** argv) {
    ArgParser a(argc,argv);

    a.addArgument("InputFileName",ArgParser::required, "Name of the Input File");
    a.addArgument("OutputFileName",ArgParser::required,"Name of the Output File");
    a.addLongOption("lumi",ArgParser::reqArg,"specify luminosity (overrides value in workspace)");

  std::string ret;
  if(a.process(ret) !=0){
    std::cout << "Invalid Options:  " << ret << std::endl;
    a.printOptions(argv[0]);
    return 0;
  }

  Fitter fits(a.getArgument("InputFileName").c_str(),a.getArgument("OutputFileName").c_str());
  if(a.longFlagPres("lumi")) fits.setLumi( atof(a.getLongFlag("lumi").c_str()) );
  fits.Run();
}
