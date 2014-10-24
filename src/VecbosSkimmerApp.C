
#include <string>
#include "TChain.h"
#include <fstream>
#include <cstring>

#include "VecbosSkimmer.hh"
#include "ArgParser.hh"
#include "PassTrigger.hh"

int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("Input",ArgParser::required,"Set the input, this is either a list of vecbos nTuples or a text file with trigger decisions");
  a.addArgument("ListOfEvents",ArgParser::required,"list of events to skim (line separated run:lumi:event)");
  a.addArgument("OutputFile",ArgParser::required,"path to the output file"); 

    string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }
  
  TChain *theChain = new TChain("ntp1");
  std::cout << "chaining" << std::endl;
  char Buffer[2000];
  TString RootFileName;
  char MyRootFile[2000];  
  ifstream *inputFile = new ifstream(a.getArgument("Input").c_str());
  if(!inputFile){
    std::cout << "Invalid input file" <<std::endl;
    return 1;
  }
  // get the tree with the conditions from the first file
  //  TTree *treeCond = new TTree();
  //  int nfiles=1;
  
  bool useXRD = a.longFlagPres("XRootD");
  TString xrdSite="";
  if(useXRD) xrdSite = a.getLongFlag("XRootD");
  
  char tmpFileName[2000];
  
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,2000);
    RootFileName = TString(Buffer).Strip();
    
    if(RootFileName.First('#') != -1 || RootFileName.Length()==0) continue;
    
    if(useXRD) {
      RootFileName.Replace(0,RootFileName.Index("/store"),
			   Form("root://xrootd.unl.edu//store/test/xrootd/%s/",xrdSite.Data()));
    }
    theChain->Add(RootFileName);
    
    std::cout << "chaining " << RootFileName << std::endl;
  } 
  inputFile->close();
  


  std::cout << "running" << std::endl;
  VecbosSkimmer sk(theChain);
  
  fstream evtFile(a.getArgument("ListOfEvents"));
  while( !(evtFile.eof()) ) {
    evtFile.getline(Buffer,2000);
    std::cout << Buffer << std::endl;
    //tokenize the string
    size_t run,lumi,evt;
    char * pch = strtok(Buffer,":");
    if(pch==NULL) continue; //not a valid line
    run = atoi(pch);
    pch = strtok(NULL,":");
    if(pch==NULL) continue; //not a valid line
    lumi = atoi(pch);
    pch = strtok(NULL,":");
    if(pch!=NULL) {
      evt = atoi(pch);
    }else{ //interpret this case as the LumiSection not being given
      evt = lumi;
      lumi = evtInfo::NOLUMI; //the default case
    }

    std::cout << run << " " << lumi << " " << evt << std::endl;
    sk.addEvent(run,lumi,evt);
  }


  sk.Skim(a.getArgument("OutputFile").c_str());
}
