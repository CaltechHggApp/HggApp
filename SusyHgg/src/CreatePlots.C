#include <iostream>
#include "../include/ArgParser.hh"
#include "Plotter.hpp"

int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("InputFile",ArgParser::required,"Name of the Input File");
  a.addArgument("OutputFolder",ArgParser::required,"Name of the folder to which to save the output");
  a.addArgument("OutputTag",ArgParser::required,"Tag to append to output plots");

  a.addLongOption("isSMS",ArgParser::noArg,"file is an SMS file (no data)");
  
  std::string ret;
  if(a.process(ret) !=0){
    std::cout << "Invalid Options:  " << ret << std::endl;
    a.printOptions(argv[0]);
    return 0;
  }

  Plotter plots(a.getArgument("InputFile"),a.getArgument("OutputFolder"),a.getArgument("OutputTag"));
  if(a.longFlagPres("isSMS")) plots.setIsSMS(true);

  plots.Run();

};
