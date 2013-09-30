#include "ArgParser.hh"
#include "MakeBiasStudy.h"
#include "ReadConfig.hh"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  ArgParser a(argc,argv);
  a.addArgument("WorkspaceFile",ArgParser::required,"path to the workspace file");
  a.addArgument("outputFile",ArgParser::required,"path to the output file");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  string wsFile = a.getArgument("WorkspaceFile");
  string outFile = a.getArgument("outputFile");
  MakeBiasStudy msht(wsFile,outFile);
  msht.run();
  return 0;

}
