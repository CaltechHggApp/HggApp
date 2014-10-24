
#include <string>
#include "TChain.h"
#include <fstream>

#include "PassTrigger.hh"
#include "ArgParser.hh"

int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("Input",ArgParser::required,"Set the input, this is either a list of vecbos nTuples or a text file with trigger decisions");
  a.addArgument("AnalysisNTuple",ArgParser::required,"analysis NTuple to check");
  a.addArgument("AnalysisNTupleTreeName",ArgParser::required,"Tree Name in the analysis NTuple to check");  
  a.addArgument("OutputFile",ArgParser::required,"path to the output file"); 

  a.addLongOption("DoOutputStep",ArgParser::noArg,"Do the output step (run on the file with trigger decisions given in Input)");


    string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }
  
  TChain *theChain = new TChain("ntp1");
  if(!(a.longFlagPres("DoOutputStep"))) {
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
  }

  TFile * analysisTreeFile = new TFile(a.getArgument("AnalysisNTuple").c_str());
  TTree * analysisTree     = (TTree*)analysisTreeFile->Get(a.getArgument("AnalysisNTupleTreeName").c_str());

  if(analysisTree==0) {
    std::cout << "ERROR: Invalid Analysis NTuple file" << std::endl;
    return -1;
  }

  std::cout << "running" << std::endl;
  PassTrigger pt(theChain,analysisTree);
  pt.addTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
  //pt.addTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v");
  //pt.addTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v");
  if( a.longFlagPres("DoOutputStep")) {
    pt.setDoOutputStep( a.getArgument("Input"));
  }

  pt.Loop(a.getArgument("OutputFile"));
  analysisTreeFile->Close();
}
