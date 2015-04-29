#include "SignalInjectionNew.hh"

SignalInjection::SignalInjection(TString outputFileName):rng(0) {
  outputFile = new TFile(outputFileName,"RECREATE");
}

SignalInjection::~SignalInjection() {
  for(auto file : smHiggsFiles) file->Close();
  dataFile->Close();
  smsFile->Close();

  outputFile->Close();
}

TH2F* SignalInjection::getToy(const TH2F& input,TString name) {

  const int nX = input.GetNbinsX();
  const int nY = input.GetNbinsY();

  if(name=="") name = input.GetName();

  TH2F* output = (TH2F*)input.Clone(name);

  for(int iX=1; iX<=nX; iX++) {
    for(int iY=0; iY<=nY; iY++) {
      output->SetBinContent(iX,iY, rng.Poisson(output->GetBinContent(iX,iY)) );
    }
  }
  return output;
}
#include <iostream>
void SignalInjection::addHiggsFile(TString name) {
  smHiggsFiles.push_back(  openInjFile(name,0) );
}

TFile* SignalInjection::openInjFile(TString fileName, TFile *file) {
  std::cout << "opening " << fileName << std::endl;
  if(file) file->Close();
  return new TFile(fileName);

  file = new TFile(fileName);

  assert(file!=0 && fileName.Data());
  std::cout << file << std::endl;
}


void SignalInjection::make() {
  assert(smsFile!=0 && smsFile->IsOpen());
  assert(dataFile!=0 && dataFile->IsOpen());
  for(auto higgsFile: smHiggsFiles) assert(higgsFile!=0 && higgsFile->IsOpen());

  for(auto cat: *cats) {
    TH2F* input_data = (TH2F*)dataFile->Get("data_"+cat+"_SidebandRegion");
    TH2F* data = getToy(*input_data,"data_"+cat+"_SignalRegion");
    for(auto higgsFile: smHiggsFiles) {
      data->Add( getToy( *(TH2F*)higgsFile->Get("data_"+cat+"_SignalRegion") ) );      
    }
    TH2F* sms_data = (TH2F*)smsFile->Get("data_"+massPoint+"_"+cat+"_SignalRegion");
    assert(xsec>=0 && xsec < 1e6);
    sms_data->Scale(xsec);
    data->Add(sms_data);
    
    outputFile->cd();
    data->Write();
  }
  TList *list = dataFile->GetListOfKeys();
  for(int i=0;i<list->GetEntries();i++) {
    TObject * o = dataFile->Get(list->At(i)->GetName());

    if( TString(o->GetName()).Contains("SignalRegion") ) continue;
    outputFile->cd();
    o->Write();
  }
}
