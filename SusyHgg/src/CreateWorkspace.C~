
#include <iostream>
#include <string>
#include "exception"
#include <stdexcept>
#include <tuple>
#include "assert.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooThresholdCategory.h"
#include "RooCategory.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TH2F.h"

#include "../include/ArgParser.hh"



std::tuple<std::string,int,std::string> parseName(std::string argument); //parses the input argument into pair(name,path), throws an exception if format is wrong

std::vector<std::pair<int,int>> getSMSPoint(std::string name); //get a list of the SMS point

void addToWorkspace(std::string name,int num_orig, std::string path, RooWorkspace* const workspace);

int main(int argc, char** argv) {
  ArgParser a(argc,argv);

  a.addArgument("OutputFileName","Name of the Output File");
  a.addUnlimitedArgument("InputFiles","list of input files (format is 'name:original_number_of_events:path_to_file')");
  a.addShortOption('a',ArgParser::noArg,"Append to existing workspace (don't recreate)");

  std::string ret;
  if(a.process(ret) !=0){
    std::cout << "Invalid Options:  " << ret << std::endl;
    a.printOptions(argv[0]);
    return 0;
  }
  
  std::string file_mode = (a.shortFlagPres('a')? "UPDATE" : "RECREATE");
  
  TFile outputFile(a.getArgument("OutputFileName").c_str(),file_mode.c_str());
  RooWorkspace *ws =0;
  TObjArray *mc_name_array = 0;

  if(a.shortFlagPres('a')) {
    ws = (RooWorkspace*)outputFile.Get("susy_hgg_workspace");
    mc_name_array = (TObjArray*)outputFile.Get("mc_name_array");
  }else{
    ws = new RooWorkspace("susy_hgg_workspace","");
    mc_name_array = new TObjArray();

  }
  if(ws==0) throw std::runtime_error("Fatal Error opening output file");

  RooCategory mc_names("mc_names","");

  for(auto input = a.getUnlimitedArgument("InputFiles"); input != a.args_end(); input++) {
    std::tuple<std::string,int,std::string> fPair = parseName(*input);
    addToWorkspace(std::get<0>(fPair),std::get<1>(fPair),std::get<2>(fPair),ws);

    if(std::get<0>(fPair).compare("data")!=0) {
      if(std::get<0>(fPair).find("sms") == std::string::npos) { //not an SMS
	mc_name_array->Add( new TObjString(std::get<0>(fPair).c_str()) );
      }else{
	std::vector<std::pair<int,int>> sms_points = getSMSPoint("");
	for(auto& pt: sms_points) {
	  mc_name_array->Add( new TObjString(Form("%s_%d_%d",std::get<0>(fPair).c_str(),pt.first,pt.second)) );
	}
      }
    }
  }
  //ws->import(mc_names);
  outputFile.cd();
  ws->Print();
  ws->Write();
  std::cout << "DONE WRITING WORKSPACE" << std::endl;
  mc_name_array->Write("mc_name_array",TObject::kSingleKey);
  outputFile.Close();
}


std::tuple<std::string,int,std::string> parseName(std::string argument) {
  //returns name,path to file
  //throws an exception if the format is wrong, but does not check that the file path is valid

  std::size_t sep = argument.find(":");
  if(sep == std::string::npos) throw std::runtime_error(Form("Invalid Argument: %s",argument.c_str()));

  std::string name = argument.substr(0,sep);

  std::string rem = argument.substr(sep+1);

  std::size_t sep2 = rem.find(":");
  if(sep2 == std::string::npos) throw std::runtime_error(Form("Invalid Argument: %s",argument.c_str()));

  std::string num = rem.substr(0,sep2);
  std::string path = rem.substr(sep2+1);

  return std::make_tuple(name,atoi(num.c_str()),path);
}


void addToWorkspace(std::string name,int num_orig, std::string path, RooWorkspace* const workspace) {
  TFile inputFile(path.c_str());
  TTree * tree= (TTree*)inputFile.Get("SusyHggTree");

  

  RooArgSet vars;
  vars.add(*(new RooRealVar("mgg","M_{#gamma#gamma}",110,150)));
  vars.add(*(new RooRealVar("ptgg","pT_{#gamma#gamma}",0,5000)));

  vars.add(*(new RooRealVar("pho1_pt","p_{T}^{#gamma-lead}"   ,0,5000)));
  vars.add(*(new RooRealVar("pho2_pt","p_{T}^{#gamma-sublead}",0,5000)));
 
  vars.add(*(new RooRealVar("pho1_eta","#eta^{#gamma-lead}",-3,3)));
  vars.add(*(new RooRealVar("pho2_eta","#eta^{#gamma-sublead}",-3,3)));

   
  vars.add(*(new RooRealVar("pho1_r9","R_{9}^{#gamma-lead}",0,1.2)));
  vars.add(*(new RooRealVar("pho2_r9","R_{9}^{#gamma-sublead}",0,1.2)));
 
  vars.add(*(new RooRealVar("pho1_sigEoE","#sigma_{E}/E^{#gamma-lead}",0,0.5)));
  vars.add(*(new RooRealVar("pho2_sigEoE","#sigma_{E}/E^{#gamma-sublead}",0,0.5)));

  vars.add(*(new RooRealVar("pho1_pass_iso","WP95 Iso lead",-0.1,1.1)));
  vars.add(*(new RooRealVar("pho2_pass_iso","WP95 Iso sublead",-0.1,1.1)));
  
  vars.add(*(new RooRealVar("ele1_pt","p_{T}^{e-lead}"  ,0,5000)));
  vars.add(*(new RooRealVar("mu1_pt" ,"p_{T}^{#mu-lead}",0,5000)));

  vars.add(*(new RooRealVar("MR","M_{R}",0,3000)));
  vars.add(*(new RooRealVar("Rsq","R^{2}",0,1)));

  vars.add(*(new RooRealVar("pileupWeight","",0,20)));

  vars.add(*(new RooRealVar("m22","",0,1000)));
  vars.add(*(new RooRealVar("m23","",0,1000)));
  vars.add(*(new RooRealVar("m24","",0,1000)));
  vars.add(*(new RooRealVar("m25","",0,1000)));

  TString puWeightString = "pileupWeight";
  if(name.compare("data")==0) puWeightString="";

  RooDataSet data(name.c_str(),"",tree,vars,"",puWeightString);

  RooRealVar *mgg = (RooRealVar*)vars.find("mgg");
  //mgg->setRange("sideband_low",120,122);
  //mgg->setRange("sideband_high",128,130);
  //mgg->setRange("signal",123,127);

  RooThresholdCategory regionNEBEB("regionNEBEB","",*(RooRealVar*)vars.find("mgg"),"Background");
  regionNEBEB.addThreshold(132,"Sideband");
  regionNEBEB.addThreshold(130,"Buffer");
  regionNEBEB.addThreshold(128,"Signal");
  regionNEBEB.addThreshold(122,"Buffer");
  regionNEBEB.addThreshold(120,"Sideband");
  regionNEBEB.addThreshold(118,"Background");

  RooThresholdCategory regionEBEB("regionEBEB","",*(RooRealVar*)vars.find("mgg"),"Background");
  regionEBEB.addThreshold(130,"Sideband");
  regionEBEB.addThreshold(128,"Buffer");
  regionEBEB.addThreshold(126.5,"Signal");
  regionEBEB.addThreshold(123.5,"Buffer");
  regionEBEB.addThreshold(122,"Sideband");
  regionEBEB.addThreshold(120,"Background");

  //data.addColumn(regionNEBEB);
  //data.addColumn(regionEBEB);


  if(name.find("sms") == std::string::npos) {
    RooRealVar N((std::string("N_")+name).c_str(),"",num_orig);
    workspace->import(data);
    workspace->import(N);
  }else{
    TFile *norm_file=0;
    if(name.find("ChiZH") !=std::string::npos) norm_file = new TFile("SMS-TChiZH_ZincHgg_2J_mChargino-130to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola__Summer12-START53_V19_FSIM-v1_ENTRIES.root");
    if(name.find("ChiWH") !=std::string::npos) norm_file = new TFile("SMS-TChiWH_WincHgg_2J_mChargino-130to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola__Summer12-START53_V19_FSIM-v1_ENTRIES.root");

    TH2F* norm_hist = 0;
    if(norm_file) norm_hist = (TH2F*)norm_file->Get("Ntotal");

    std::vector<std::pair<int,int>> sms_points = getSMSPoint("");
    for(auto& pt: sms_points) {
      int Num = (norm_hist==0 ? 1 : norm_hist->GetBinContent(norm_hist->FindFixBin(pt.second,pt.first)));
      RooRealVar N(Form("N_%s_%d_%d",name.c_str(),pt.first,pt.second),"",Num);

      RooDataSet * data_pt = (RooDataSet*)data.reduce(Form("abs(m22-%d)<10 && abs(m23-%d)<10",pt.first,pt.second));
      workspace->import(*data_pt,RooFit::Rename( Form("%s_%d_%d",name.c_str(),pt.first,pt.second)) );
      workspace->import(N);
    }
    delete norm_file;
  }
  
}


std::vector<std::pair<int,int>> getSMSPoint(std::string name) {
  std::vector<std::pair<int,int>> res;
  for(int iM23=125;iM23<=500;iM23+=25) {
    for(int iM22=0;iM22<=iM23-125;iM22+=25) {
      res.push_back(std::make_pair(iM22,iM23));
    }
  }
  return res;
}
