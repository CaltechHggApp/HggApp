
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

    a.addArgument("InputFileName", "Name of the Input File");
    a.addArgument("OutputFileName","Name of the Output File");

  std::string ret;
  if(a.process(ret) !=0){
    std::cout << "Invalid Options:  " << ret << std::endl;
    a.printOptions(argv[0]);
    return 0;
  }

  Fitter fits(a.getArgument("InputFileName").c_str(),a.getArgument("OutputFileName").c_str());

  fits.Run();
}
