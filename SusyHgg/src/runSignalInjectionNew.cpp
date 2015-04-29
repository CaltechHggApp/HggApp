/**

creates a signal injection

*/


#include <iostream>
#include <string>

#include "ArgParser.hh"

#include "SignalInjectionNew.hh"

using namespace std;

int main(int argc,char** argv) {
  ArgParser a(argc, argv);

  a.addArgument("Folder",ArgParser::required,"folder with the data and to which to write the signal injection");
  a.addArgument("SMSName",ArgParser::required,"type of SMS to inject (same as file name)");
  a.addArgument("massPoint",ArgParser::required,"mass point to inject");
  a.addLongOption("xsec",ArgParser::reqArg,"x-sec to inject");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }
  
  string folder = a.getArgument("Folder");
  cout << folder << endl;
  string smsName = a.getArgument("SMSName");
  string massPoint = a.getArgument("massPoint");

  float xsec =1;
  if(a.longFlagPres("xsec")) xsec = atof(a.getLongFlag("xsec").c_str());
  SignalInjection injector(folder+"/signalInj_"+smsName+"_"+massPoint+"_"+std::to_string(xsec)+"pb.root");
  //cout << injector.getDataFile() << endl;
  injector.setDataFile(folder+"/data.root");
  //cout << injector.getDataFile() << endl;
  injector.addHiggsFile(folder+"/ggH.root");
  injector.addHiggsFile(folder+"/vbfH.root");
  injector.addHiggsFile(folder+"/wzH.root");
  injector.setSMSFile(folder+"/"+smsName+".root");

  injector.setMassPoint(massPoint);
  injector.setInjXsec( xsec );

  injector.make();

}
