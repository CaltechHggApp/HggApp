#include "ArgParser.hh"
#include "MakeSpinFits.h"
#include "ReadConfig.hh"

#include <iostream>
#include <string>
#include <map>

int main(int argc, char** argv){
  using namespace std;
  
  ArgParser a(argc,argv);
  a.addArgument("InputWorkspace",ArgParser::required,"input workspace");
  a.addArgument("OutputWorkspace",ArgParser::required,"output workspace");
  a.addLongOption("BkgFit",ArgParser::reqArg,"Background Fit Type [poly,sexp,dexp,texp,pow,dpow] (default: dexp)");
  a.addLongOption("MCSamples",ArgParser::reqArg,"Specify the MC sample to process (comma separated) or none (default: all)");
  a.addLongOption("CrystalBall",ArgParser::noArg,"convolute a crystal ball with the signal model");
  a.addLongOption("DY",ArgParser::noArg,"Specify that we are fitting the Zee peak");
  a.addLongOption("FixPeak",ArgParser::noArg,"fix the peak position of the signal model at the generated mass");
  a.addLongOption("BkgOnly",ArgParser::noArg,"only run the background fits");

  string ret;
  if(a.process(ret) != 0){
    a.printOptions(argv[0]);
    return 0;
  }

  string inputWS = a.getArgument("InputWorkspace");
  string outputWS = a.getArgument("OutputWorkspace");
  string fit = "exp";
  if(a.longFlagPres("BkgFit")){
    fit = a.getLongFlag("BkgFit");    
  }

  std::map<string,MakeSpinFits::BkgFitType> validFits = {
      { "poly", MakeSpinFits::kPoly },
      { "pow",  MakeSpinFits::kPow },
      { "dpow", MakeSpinFits::kDoublePow },
      { "sexp", MakeSpinFits::kSingleExp },
      { "dexp", MakeSpinFits::kDoubleExp },
      { "texp", MakeSpinFits::kTripleExp },
      { "mexp", MakeSpinFits::kModifiedExp }
  };
  
  MakeSpinFits msf(inputWS,outputWS);
  auto matchedFit = validFits.find(fit);
  if( matchedFit == validFits.end()) {
      std::cout << "Invalid Fit " << fit << std::endl;
      return -1;
  }
  msf.setBkgFit(matchedFit->second);
  if(a.longFlagPres("MCSamples")) // specify the samples to test
    {
      vector<string> samples = ReadConfig::tokenizeString( a.getLongFlag("MCSamples"), ",");
      if(samples.at(0).compare("none")==0) samples.clear();
      if(msf.specifySamples(samples) != 0) 
	{//error in a sample specification
	  return -1;
	}
    }

  if(a.longFlagPres("CrystalBall")) msf.setUseCrystalBall();
  if(a.longFlagPres("DY")) msf.setMeanRange(88,94,91);

  if(a.longFlagPres("FixPeak")){
    if(a.longFlagPres("DY")) msf.setMeanRange(90,92,91);
    else msf.setMeanRange(124,126,125);
  }

  if(a.longFlagPres("BkgOnly")) msf.setMakeBkgOnly();

  msf.setAddSWeight(true);
  msf.run();
  msf.save();

}
