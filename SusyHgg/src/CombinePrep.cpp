#include "CombinePrep.hh"
#include <sys/stat.h>

#include "FitterNew.hpp"
#include "SMSFitterNew.hh"

void CombinePrep::openFiles() {
  dataFile = std::unique_ptr<TFile>(new TFile(dataFilePath));
  assert(dataFile->IsOpen() && Form("could not open data file %s",dataFilePath.Data()));

  for(auto smHiggs: smHiggsFilePaths) {
    smHiggsFiles[smHiggs.first] = std::unique_ptr<TFile>(new TFile(smHiggs.second));
    assert(dataFile->IsOpen() && Form("could not open sm higgs file %s",smHiggs.second.Data()));    
  }

  for(auto sms: smsFilePaths) {
    smsFiles[sms.first] = std::unique_ptr<TFile>(new TFile(sms.second));
    assert(dataFile->IsOpen() && Form("could not open sms file %s",sms.second.Data()));    
  }
}

void CombinePrep::Make(){
  makeRootFile();
  std::vector<TString> sms_points;
  SMSFitter::getSMSPoints(&sms_points);
  for(auto pt: sms_points) {
    makeOnePoint("sms_ChiWH",pt);
    makeOnePoint("sms_ChiZH",pt);
  }
}

void CombinePrep::makeOnePoint(TString smsName, TString smsPoint) {
  std::fstream dataCardFile( Form("%s/combine/%s_%s.txt",outputFolder.Data(),smsName.Data(),smsPoint.Data()), std::fstream::out);

  dataCardFile << "imax\t1\njmax\t*\nkmax\t*\n";
  dataCardFile << "-------------------------------------\n";
  dataCardFile << "bin 1\n";
  
  dataCardFile << "observation\t" << yields["data"] << std::endl;
  dataCardFile << "-------------------------------------\n";
  dataCardFile << "shapes * * data.root $PROCESS $PROCESS_$SYSTEMATIC\n";
  dataCardFile << "-------------------------------------\n";
  
  dataCardFile << "bin\t\t\t1\t\t1"; //signal     background     
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t1";
  dataCardFile << "\n";

  dataCardFile << "process\t\t\t"<< smsName+"_"+smsPoint <<"\tbkg";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    dataCardFile << "\t\t" << sm->first;
  }
  dataCardFile << "\n";

  dataCardFile << "process\t\t\t0\t\t1"; //signal     background     
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t" << i+2;
  dataCardFile << "\n";
  
  dataCardFile << "rate\t\t\t" << yields[smsName+"_"+smsPoint] << "\t\t" << yields["bkg"] ;
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    dataCardFile << "\t\t" << yields[sm->first];
  }
  dataCardFile << "\n";

  dataCardFile << "-------------------------------------\n";

  //systematics!!

  //lumi
  dataCardFile << "lumi\tlnN\t\t\t1.025\t\t-";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t1.025";
  dataCardFile << "\n";

  //Higgs scales
  
  dataCardFile << "QCDScale_VH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "wzH") dataCardFile << "\t\t0.982/1.021";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "QCDScale_qqH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "vbfH") dataCardFile << "\t\t0.992/1.003";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "QCDScale_ggH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "ggH") dataCardFile << "\t\t0.918/1.076";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "pdf_VH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "wzH") dataCardFile << "\t\t0.958/1.042";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "pdf_qqH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "vbfH") dataCardFile << "\t\t0.972/1.026";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "pdf_ggH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "ggH") dataCardFile << "\t\t0.930/1.076";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";


  //fit systematic
  dataCardFile << "fit\tshape\t\t\t-\t\t1";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
  dataCardFile << "\n";
    
  //bkgShape systematic
  dataCardFile << "bkgShape\tshape\t\t\t-\t\t1";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
  dataCardFile << "\n";
    
  //mc systematics
  for(auto sys: *sysNames) {
    dataCardFile << sys << "\tshape\t\t\t1\t\t-";
    for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t1";
    dataCardFile << "\n";
  }


  dataCardFile.close();

}

void CombinePrep::makeRootFile() {
  int dir = mkdir( Form("%s/combine/",outputFolder.Data()), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  assert((dir==0 || errno == EEXIST) && Form("could not create combine directory %s/combine/",outputFolder.Data()));
  if(dataFile == nullptr) openFiles();

  TString outputFileStub = outputFolder+"/combine/data"; // <folder>/combine/<name>_<pt>  .txt and .root will be added as appropriate

  TFile outputRootFile(outputFileStub+".root","RECREATE");

  TH1F* obsHist = makeCategoryHistogram(dataFile.get(),"data_obs","SignalRegion");

  yields["data"] = obsHist->Integral();

  TH1F* bkgHist = makeCategoryHistogram(dataFile.get(),"bkg","SidebandRegion");

  yields["bkg"] = bkgHist->Integral();
  TH1F* bkgUpHist = makeCategoryHistogram(dataFile.get(),"bkg_fitUp","SidebandRegion_fitUp");
  TH1F* bkgDownHist = makeCategoryHistogram(dataFile.get(),"bkg_fitDown","SidebandRegion_fitDown");
  makeCategoryHistogram(dataFile.get(),"bkg_bkgShapeUp","SidebandRegion_bkgShapeUp")->Write();
  makeCategoryHistogram(dataFile.get(),"bkg_bkgShapeDown","SidebandRegion_bkgShapeDown")->Write();
  

  std::cout << "Adding SM Higgs" << std::endl;
  for(auto sm = smHiggsFiles.begin(); sm !=smHiggsFiles.end(); sm++) {
    TH1F* h = makeCategoryHistogram(sm->second.get(),sm->first,"SignalRegion");
    yields[sm->first] = h->Integral();
    h->Write();
    std::cout << sm->first << std::endl;
    for(auto sys: *sysNames) {
      std::cout << sys << std::endl;
      makeCategoryHistogram(sm->second.get(),sm->first+"_"+sys+"Up","SignalRegion_"+sys+"_Up")->Write();     
      makeCategoryHistogram(sm->second.get(),sm->first+"_"+sys+"Down","SignalRegion_"+sys+"_Down")->Write();     
    }
  }

  std::vector<TString> sms_points;
  SMSFitter::getSMSPoints(&sms_points);
  
  for(auto sms = smsFiles.begin(); sms !=smsFiles.end(); sms++) {
    for(auto pt : sms_points) {
      TH1F* h = makeCategoryHistogram(sms->second.get(),sms->first+"_"+pt,"SignalRegion",pt);
      yields[sms->first+"_"+pt] = h->Integral();
      h->Write();
      for(auto sys: *sysNames) {
	makeCategoryHistogram(sms->second.get(),sms->first+"_"+pt+"_"+sys+"Up","SignalRegion_"+sys+"_Up",pt)->Write();     
	makeCategoryHistogram(sms->second.get(),sms->first+"_"+pt+"_"+sys+"Down","SignalRegion_"+sys+"_Down",pt)->Write();     
      }
    }
  }
  
  



  obsHist->Write();
  bkgHist->Write();
  bkgUpHist->Write();
  bkgDownHist->Write();

  outputRootFile.Close();
}

std::unique_ptr<TH1F> CombinePrep::compressHistogram(const TH2F& input) {
  int nX = input.GetNbinsX();
  int nY = input.GetNbinsY();
  TString name = input.GetName();

  std::unique_ptr<TH1F> compressed(new TH1F(name,"",nX*nY,-0.5,nX*nY-0.5));
  
  for(int iX=1;iX<=nX;iX++) {
    for(int iY=1;iY<=nY;iY++) {
      compressed->SetBinContent((iX-1)*(nY+1)+iY,input.GetBinContent(iX,iY));
    }
  }
  return compressed;
}

TH1F* CombinePrep::makeCategoryHistogram(TFile *file,TString histName,TString postfix,TString sms_pt) {
  
  if(sms_pt!="")sms_pt+="_";

  //std::cout << file->GetName() << std::endl;


  TH2F* dummy = (TH2F*)file->Get( "data_"+sms_pt+*(catNames->begin())+"_"+postfix );

  assert(dummy != 0 && "could not get histogram from file to make categories");

  int nX=dummy->GetNbinsX();
  int nY=dummy->GetNbinsY();

  TH1F* catHist = new TH1F(histName,"",nX*nY*catNames->size(),-0.5,nX*nY*catNames->size()-0.5);

  for(int iCat=0;iCat<catNames->size();iCat++) {
    TH2F* hist = (TH2F*)file->Get( "data_"+sms_pt+catNames->at(iCat)+"_"+postfix );
    std::unique_ptr<TH1F> compHist = compressHistogram(*hist);
    //std::cout << histName << " " << postfix << " " << sms_pt << " " << iCat << "    " << dummy->Integral() << " --- " <<  compHist->Integral() << std::endl;
    for(int iBin=1;iBin<=nX*nY;iBin++) {
      catHist->SetBinContent(iBin+nX*nY*iCat,compHist->GetBinContent(iBin));
    }
  }
  return catHist;
}
