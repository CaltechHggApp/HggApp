#include "CombinePrep.hh"
#include <sys/stat.h>

#include "FitterNew.hpp"
#include "SMSFitterNew.hh"
#include <algorithm>

#define UNCORR_DIFF 1

#define MAX(a,b) (((a)>(b))?(a):(b))

void CombinePrep::openFiles() {
  dataFile = std::unique_ptr<TFile>(new TFile(dataFilePath));
  std::cout << "data file: " << dataFilePath << std::endl;
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
  TString outputSubDir = dataFilePath;
  outputSubDir.Remove(0,outputSubDir.Last('/')+1);
  outputSubDir.Remove(outputSubDir.Last('.'));
  outputFolder = outputFolder +"/"+outputSubDir+"/";

  makeRootFile();

  //hist_nBins = nX*nY*catNames->size();
  if(useVarBinning) {
    hist_nBins=0;
    for(auto binning = binningMap.begin(); binning != binningMap.end();binning++) hist_nBins+= (binning->second.size()-1)*binning->second.size()/2;
  }

  std::vector<TString> sms_points;
  SMSFitter::getSMSPoints(&sms_points,isFullSIM);
  for(auto pt: sms_points) {
    makeOnePoint("sms_ChiWH",pt);
    //makeOnePoint("sms_ChiZH",pt);
    //makeOnePoint("sms_ChiWH_FSIMSmear",pt);
    //makeOnePoint("sms_ChiZH_FSIMSmear",pt);
    makeOnePoint("sms_ChiHH",pt);
    //makeOnePoint("sms_ChiHH_FSIMSmear",pt);
  }
}

void CombinePrep::makeOnePoint(TString smsName, TString smsPoint) {
  std::fstream dataCardFile( Form("%s/%s_%s.txt",outputFolder.Data(),smsName.Data(),smsPoint.Data()), std::fstream::out);

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

  //trigger
  dataCardFile << "trigger\tlnN\t\t\t"<< triggerEffErr+1 <<"\t\t-";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t" << triggerEffErr+1;
  dataCardFile << "\n";
  

  //Higgs Norm Uncertainty --NOT USED
  //dataCardFile << "HiggsScale\tlnN\t\t\t-\t\t-";
  //for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
  //dataCardFile << "\t\t0.74/1.28";
  //}
  //dataCardFile << "\n";

  //Higgs scales  
  dataCardFile << "QCDScale_VH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "wzH") dataCardFile << "\t\t0.982/1.021";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "QCDScale_qqH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "vbfH") dataCardFile << "\t\t0.998/1.002";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "QCDScale_ggH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "ggH") dataCardFile << "\t\t0.928/1.072";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "QCDScale_ttH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first=="ttH") dataCardFile << "\t\t0.915/1.038";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  dataCardFile << "pdf_qqbar\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "wzH") dataCardFile << "\t\t0.977/1.024";
    else if( sm->first == "vbfH") dataCardFile << "\t\t0.973/1.026";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";

  /*
  dataCardFile << "pdf_qqH\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "vbfH") dataCardFile << "\t\t0.972/1.026";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";
  */
  dataCardFile << "pdf_gg\tlnN\t\t\t-\t\t-";
  for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
    if( sm->first == "ggH") dataCardFile << "\t\t0.935/1.075";
    else if( sm->first=="ttH") dataCardFile << "\t\t1.081/0.925";
    else dataCardFile << "\t\t-";
  }
  dataCardFile << "\n";


#if UNCORR_DIFF
  //bkgShape systematic
  //dataCardFile << "bkgShape\tshape\t\t\t-\t\t1";
  //for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
  //dataCardFile << "\n";

  //bkgstatistics systematic -- uncorrelated in each bin
  for(int iBin=0;iBin<hist_nBins; iBin++) {
    dataCardFile << "bkgShape_bin" << iBin << "\tshape\t\t\t-\t\t1";
    for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
    dataCardFile << "\n";
  }

#else
  //bkgShape systematic
  dataCardFile << "bkgShape\tshape\t\t\t-\t\t1";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
  dataCardFile << "\n";

  //scale factor systematic
  dataCardFile << "sfSyst\tshape\t\t\t-\t\t1";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
  dataCardFile << "\n";
#endif
  /*
  //bkgstatistics systematic
  dataCardFile << "sidebandStatistics\tshape\t\t\t-\t\t1";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
  dataCardFile << "\n";
  */

  //bkgstatistics systematic -- uncorrelated in each bin
  for(int iBin=0;iBin<hist_nBins; iBin++) {
    dataCardFile << "sidebandStatistics_bin" << iBin << "\tshape\t\t\t-\t\t1";
    for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
    dataCardFile << "\n";
  }
    
  //fit systematic
  dataCardFile << "fit\tshape\t\t\t-\t\t1";
  for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t-";
  dataCardFile << "\n";
    
    
  //mc systematics
  for(auto sys: *sysNames) {
    if(sys=="statistics") {
      // sms statistics
      for(int iBin=0;iBin<hist_nBins; iBin++) {
	dataCardFile << sys << "_bin" << iBin << "\tshape\t\t\t1\t\t-";
	for(int i=0;i<smHiggsFiles.size();i++) {
	  dataCardFile << "\t\t-";
	}
	dataCardFile << "\n";      
      }
      for(auto sm = smHiggsFiles.begin(); sm!= smHiggsFiles.end(); sm++) {
	for(int iBin=0;iBin<hist_nBins; iBin++) {
	  dataCardFile << sys << "_" << sm->first << "_bin" << iBin << "\tshape\t\t\t-\t\t-";
	  for(auto sm2 = smHiggsFiles.begin(); sm2!= smHiggsFiles.end(); sm2++) {
	    if(sm2->first==sm->first) dataCardFile << "\t\t1";
	    else dataCardFile << "\t\t-";
	  }
	  dataCardFile << "\n";      
	}
      }
    }else{
      dataCardFile << sys << "\tshape\t\t\t1\t\t-";
      for(int i=0;i<smHiggsFiles.size();i++) dataCardFile << "\t\t1";
      dataCardFile << "\n";
    }
  }


  dataCardFile.close();

}

void CombinePrep::makeRootFile() {
  int dir = mkdir( Form("%s/",outputFolder.Data()), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  assert((dir==0 || errno == EEXIST) && Form("could not create combine directory %s/combine/",outputFolder.Data()));
  if(dataFile == nullptr) openFiles();

  if(useVarBinning) {
    for(auto cat: *catNames) {
      std::cout << "--------------------------------\n" << cat << "\n--------------------------------\n" << std::endl;
      if(!externalBinning) {
	binningMap.insert( std::pair<TString,std::vector<int>>(cat,{}) ); //setup the binning map
	TH2F* totalBkg = (TH2F*)dataFile->Get("data_"+cat+"_SidebandRegion")->Clone("tmp");
	for(auto h = smHiggsFiles.begin(); h!=smHiggsFiles.end(); h++) {
	  totalBkg->Add( (TH2F*)h->second->Get("data_"+cat+"_SignalRegion") );
	}
	defineBinning(*totalBkg, binningMap[cat],1,1,1.);
	delete totalBkg;
      } else {
	for(int i=0;i<binningMap[cat].size();i++) std::cout << binningMap[cat].at(i) << " ";
	std::cout << std::endl;
      }
    }
  }

  ofstream binningFile( Form("%s/binningMap.txt",outputFolder.Data()) );
  for(auto cat: *catNames) {
    binningFile << cat << " ";
    for(auto b: binningMap[cat]) binningFile << b << " ";
    binningFile << std::endl;
  }
  binningFile.close();

  TString outputFileStub = outputFolder+"/data"; // <folder>/<name>_<pt>  .txt and .root will be added as appropriate

  TFile outputRootFile(outputFileStub+".root","RECREATE");

  TH1F* obsHist = makeCategoryHistogram(dataFile.get(),"data_obs","SignalRegion");

  yields["data"] = obsHist->Integral();

  TH1F* bkgHist = makeCategoryHistogram(dataFile.get(),"bkg","SidebandRegion");

  yields["bkg"] = bkgHist->Integral();
  TH1F* bkgUpHist = makeCategoryHistogram(dataFile.get(),"bkg_fitUp","SidebandRegion_fitUp");
  TH1F* bkgDownHist = makeCategoryHistogram(dataFile.get(),"bkg_fitDown","SidebandRegion_fitDown");

  TH1F* statistics = makeCategoryHistogram(dataFile.get(),"bkg_statistics","SidebandRegion_statistics");

  //compute the per-bin scale factors
  TH1F* scaleFactors = (TH1F*)bkgHist->Clone("scaleFactors");
  scaleFactors->Divide(statistics);
  
  TH1F* statisticsUp = (TH1F*)statistics->Clone("statisticsUp");
  TH1F* statisticsDown = (TH1F*)statistics->Clone("statisticsDown");

  for(int i=1;i<=statistics->GetNbinsX();i++) {
    float N = statistics->GetBinContent(i);
    statisticsUp->SetBinContent(i,(N+sqrt(N))*scaleFactors->GetBinContent(i));
    statisticsDown->SetBinContent(i,(N-sqrt(N))*scaleFactors->GetBinContent(i));
  }

  std::vector<TH1F*> sidestatUp   = uncorrelateHistogram(statisticsUp,bkgHist,"bkg_sidebandStatistics_bin%dUp"); 
  std::vector<TH1F*> sidestatDown = uncorrelateHistogram(statisticsDown,bkgHist,"bkg_sidebandStatistics_bin%dDown"); 


  //std::vector<TH1F*> sidestatUp = makeUncorrPerBinCategoryHistogram(dataFile.get(),"bkg_sidebandStatistics_bin%dUp","SidebandRegion_statisticsUp","SidebandRegion");
  //std::vector<TH1F*> sidestatDown = makeUncorrPerBinCategoryHistogram(dataFile.get(),"bkg_sidebandStatistics_bin%dDown","SidebandRegion_statisticsDown","SidebandRegion");

  scaleFactors->Write();
  statistics->Write();
  for(auto h: sidestatUp) h->Write();
  for(auto h: sidestatDown) h->Write();


#if UNCORR_DIFF
  TH1F* statistics_High = makeCategoryHistogram(dataFile.get(),"bkg_bkgShapeHigh","SidebandRegion_statistics_High");
  TH1F* statistics_Low = makeCategoryHistogram(dataFile.get(),"bkg_bkgShapeLow","SidebandRegion_statistics_Low");

  TH1F* bkgShapeUp = (TH1F*)statistics_High->Clone("bkgShapeUp");
  TH1F* bkgShapeDown = (TH1F*)statistics_Low->Clone("bkgShapeDown");

  for(int i=1;i<=statistics->GetNbinsX();i++) {
    float N = statistics->GetBinContent(i);    
    float N_Low = statistics_Low->GetBinContent(i);
    float N_High = statistics_High->GetBinContent(i);
    if(N_High<N_Low) std::swap(N_High,N_Low);

    //if the difference is small we are limited by the statistical error
    if(N_High-N_Low<sqrt(N)) {
      bkgShapeUp->SetBinContent(i,(N+sqrt(N))*scaleFactors->GetBinContent(i));
      bkgShapeDown->SetBinContent(i,(N-sqrt(N))*scaleFactors->GetBinContent(i));
    } else {
      float sf_up   = scaleFactors->GetBinContent(i)*N/N_High;
      bkgShapeUp->SetBinContent(i,N*scaleFactors->GetBinContent(i) + sqrt(N_High)*sf_up);
      if(N_Low>0) {
	float sf_down = scaleFactors->GetBinContent(i)*N/N_Low;
	bkgShapeDown->SetBinContent(i,N*scaleFactors->GetBinContent(i) + sqrt(N_Low)*sf_down);
      } else {
	bkgShapeDown->SetBinContent(i,0);
      }
    }
  }

  std::vector<TH1F*> ucBkgShapeUp   = uncorrelateHistogram(bkgShapeUp,bkgHist,"bkg_bkgShape_bin%dUp"); 
  std::vector<TH1F*> ucBkgShapeDown = uncorrelateHistogram(bkgShapeDown,bkgHist,"bkg_bkgShape_bin%dDown"); 
  //std::vector<TH1F*> bkgShapeUp = makeUncorrPerBinCategoryHistogram(dataFile.get(),"bkg_bkgShape_bin%dUp","SidebandRegion_bkgShapeUp","SidebandRegion");
  //std::vector<TH1F*> bkgShapeDown = makeUncorrPerBinCategoryHistogram(dataFile.get(),"bkg_bkgShape_bin%dDown","SidebandRegion_bkgShapeDown","SidebandRegion");
  statistics_High->Write();
  statistics_Low->Write();
  for(auto h: ucBkgShapeUp) h->Write();
  for(auto h: ucBkgShapeDown) h->Write();
#else
  makeCategoryHistogram(dataFile.get(),"bkg_bkgShapeUp","SidebandRegion_bkgShapeUp")->Write();
  makeCategoryHistogram(dataFile.get(),"bkg_bkgShapeDown","SidebandRegion_bkgShapeDown")->Write();
#endif

  //makeCategoryHistogram(dataFile.get(),"bkg_sidebandStatisticsUp","SidebandRegion_statisticsUp")->Write();
  //makeCategoryHistogram(dataFile.get(),"bkg_sidebandStatisticsDown","SidebandRegion_statisticsDown")->Write();
  
  makeCategoryHistogram(dataFile.get(),"bkg_sfSystUp","SidebandRegion_sfSystUp")->Write();
  makeCategoryHistogram(dataFile.get(),"bkg_sfSystDown","SidebandRegion_sfSystDown")->Write();
  

  std::cout << "Adding SM Higgs" << std::endl;
  for(auto sm = smHiggsFiles.begin(); sm !=smHiggsFiles.end(); sm++) {
    TH1F* h = makeCategoryHistogram(sm->second.get(),sm->first,"SignalRegion");

    yields[sm->first] = h->Integral();
    h->Write();
    std::cout << sm->first << std::endl;
    for(auto sys: *sysNames) {
      std::cout << sys << std::endl;
      if(sys=="statistics") {
	std::vector<TH1F*> Up = makeUncorrPerBinCategoryHistogram(sm->second.get(),sm->first+"_"+sys+"_"+sm->first+"_bin%dUp","SignalRegion_"+sys+"_Up","SignalRegion");
	std::vector<TH1F*> Down = makeUncorrPerBinCategoryHistogram(sm->second.get(),sm->first+"_"+sys+"_"+sm->first+"_bin%dDown","SignalRegion_"+sys+"_Down","SignalRegion");

	for(auto h: Up) h->Write();
	for(auto h: Down) h->Write();


      } else{
	makeCategoryHistogram(sm->second.get(),sm->first+"_"+sys+"Up","SignalRegion_"+sys+"_Up")->Write();     
	makeCategoryHistogram(sm->second.get(),sm->first+"_"+sys+"Down","SignalRegion_"+sys+"_Down")->Write();     
      }
    }
  }

  std::vector<TString> sms_points;
  SMSFitter::getSMSPoints(&sms_points,isFullSIM);
  
  for(auto sms = smsFiles.begin(); sms !=smsFiles.end(); sms++) {
    for(auto pt : sms_points) {
      TH1F* h = makeCategoryHistogram(sms->second.get(),sms->first+"_"+pt,"SignalRegion",pt);
      yields[sms->first+"_"+pt] = h->Integral();
      h->Write();
      for(auto sys: *sysNames) {
	if(sys=="statistics") {
	  std::vector<TH1F*> Up = makeUncorrPerBinCategoryHistogram(sms->second.get(),sms->first+"_"+pt+"_"+sys+"_bin%dUp","SignalRegion_"+sys+"_Up","SignalRegion",pt);
	  std::vector<TH1F*> Down = makeUncorrPerBinCategoryHistogram(sms->second.get(),sms->first+"_"+pt+"_"+sys+"_bin%dDown","SignalRegion_"+sys+"_Down","SignalRegion",pt);

	  for(auto h: Up) h->Write();
	  for(auto h: Down) h->Write();


	} else{
	  makeCategoryHistogram(sms->second.get(),sms->first+"_"+pt+"_"+sys+"Up","SignalRegion_"+sys+"_Up",pt)->Write();     
	  makeCategoryHistogram(sms->second.get(),sms->first+"_"+pt+"_"+sys+"Down","SignalRegion_"+sys+"_Down",pt)->Write();     
	}
      }
    }
  }
  
  



  obsHist->Write();
  bkgHist->Write();
  bkgUpHist->Write();
  bkgDownHist->Write();

  outputRootFile.Close();
}

std::unique_ptr<TH1F> CombinePrep::compressHistogram(const TH2F& input,bool varBinning,int minWidth,int minBin,std::vector<int> *xBinEdges) {
  assert( !(varBinning && xBinEdges==0) );

  int nX = input.GetNbinsX();
  int nY = input.GetNbinsY();
  TString name = input.GetName();

  int nBins = nX*nY;
  if(varBinning) nBins = (xBinEdges->size()-1)*xBinEdges->size()/2;

  std::unique_ptr<TH1F> compressed(new TH1F(name,"",nBins,-0.5,nBins-0.5));

  if(varBinning) {    
    int nXBins = xBinEdges->size()-1;
    int index=1;
    for(int iX=0;iX<nXBins;iX++) {
      for(int iY=1;iY< iX+1;iY++) {
// 	std::cout << xBinEdges->at(iX) << "--" << xBinEdges->at(iX+1) << std::endl;
// 	std::cout << iY << std::endl;
// 	std::cout << input.Integral(xBinEdges->at(iX),xBinEdges->at(iX+1),iY,iY) << std::endl;
	compressed->SetBinContent( index++, input.Integral(xBinEdges->at(iX+1),xBinEdges->at(iX)-1,iY,iY) );
      }
      compressed->SetBinContent( index++, input.Integral(xBinEdges->at(iX+1),xBinEdges->at(iX)-1,iX+1,nY+1) );
    }    

    if(index != nBins+1) {
      std::cout << index << "  " << nBins+1 << std::endl;
      assert(false);
    }
  }else{
    for(int iX=minBin;iX<=nX;iX++) {
      for(int iY=1;iY<=nY;iY++) {
	compressed->SetBinContent((iX-1)*(nY+1)+iY,input.GetBinContent(iX,iY));
      }
    }
  }
  //std::cout << "bin 1:   " << compressed->GetBinContent(1) << std::endl;
  return compressed;
}

TH1F* CombinePrep::makeCategoryHistogram(TFile *file,TString histName,TString postfix,TString sms_pt) {
  
  if(sms_pt!="")sms_pt+="_";

  //std::cout << file->GetName() << std::endl;


  TH2F* dummy = (TH2F*)file->Get( "data_"+sms_pt+*(catNames->begin())+"_"+postfix );

  std::cout << "data_"+sms_pt+*(catNames->begin())+"_"+postfix << std::endl;
  assert(dummy != 0 && "could not get histogram from file to make categories");

  int nX=dummy->GetNbinsX();
  int nY=dummy->GetNbinsY();

  int nBins = nX*nY*catNames->size();
  if(useVarBinning) {
    nBins=0;
    for(auto binning = binningMap.begin(); binning != binningMap.end();binning++) nBins+= (binning->second.size()-1)*binning->second.size()/2;
  }

  TH1F* catHist = new TH1F(histName,"",nBins,-0.5,nX*nY*catNames->size()-0.5);

  int offset=0;
  for(int iCat=0;iCat<catNames->size();iCat++) {
    //std::cout << iCat << "   " << offset << std::endl;
    TH2F* hist = (TH2F*)file->Get( "data_"+sms_pt+catNames->at(iCat)+"_"+postfix );
    std::unique_ptr<TH1F> compHist = compressHistogram(*hist,useVarBinning,1,1,&binningMap[catNames->at(iCat)]);
    //std::cout << histName << " " << postfix << " " << sms_pt << " " << iCat << "    " << dummy->Integral() << " --- " <<  compHist->Integral() << std::endl;
    int nBinsThisCat=nX*nY;
    if(useVarBinning) {      
      int nXBins = binningMap[ catNames->at(iCat) ].size();
      nBinsThisCat=nXBins*(nXBins-1)/2;
    }
    for(int iBin=1;iBin<=nBinsThisCat;iBin++) {
      catHist->SetBinContent(iBin+offset,compHist->GetBinContent(iBin));
    }
    offset+=nBinsThisCat;
  }
  //std::cout << ">>>>>bin 1:   " << catHist->GetBinContent(1) << std::endl;
  return catHist;
}

std::vector<TH1F*> CombinePrep::makeUncorrPerBinCategoryHistogram(TFile *file,TString histName,TString postfix,TString nom_postfix, TString sms_pt) {

  TH1F* correlated = makeCategoryHistogram(file,"TMP",postfix,sms_pt);
  TH1F* nominal    = makeCategoryHistogram(file,"TMPNOM",nom_postfix,sms_pt);

  std::vector<TH1F*> uncorr =  uncorrelateHistogram(correlated,nominal,histName);

  delete correlated;
  return uncorr;
  
  if(sms_pt!="")sms_pt+="_";

  //std::cout << file->GetName() << std::endl;


  TH2F* dummy = (TH2F*)file->Get( "data_"+sms_pt+*(catNames->begin())+"_"+postfix );

  std::cout << "data_"+sms_pt+*(catNames->begin())+"_"+postfix << std::endl;
  assert(dummy != 0 && "could not get histogram from file to make categories");

  int nX=dummy->GetNbinsX();
  int nY=dummy->GetNbinsY();

  int nBins = nX*nY*catNames->size();
  if(useVarBinning) {
    nBins=0;
    for(auto binning = binningMap.begin(); binning != binningMap.end();binning++) nBins+= (binning->second.size()-1)*binning->second.size()/2;
  }


  std::vector<TH1F*> catHists(nBins);
  for(int i=0; i<nBins; i++) {
    catHists[i] = new TH1F(Form(histName.Data(),i),"",nBins,-0.5,nX*nY*catNames->size()-0.5);
  }
  //first, make it each histogram a copy of the nominal
  int offset=0;
  for(int iCat=0;iCat<catNames->size();iCat++) {
    //std::cout << iCat << "   " << offset << std::endl;
    TH2F* nominal = (TH2F*)file->Get( "data_"+sms_pt+catNames->at(iCat)+"_"+nom_postfix );
    std::unique_ptr<TH1F> compHist = compressHistogram(*nominal,useVarBinning,1,1,&binningMap[catNames->at(iCat)]);
    int nBinsThisCat=nX*nY;
    if(useVarBinning) {      
      int nXBins = binningMap[ catNames->at(iCat) ].size();
      nBinsThisCat=nXBins*(nXBins-1)/2;
    }
    for(int iHist=0; iHist<nBins; iHist++) {
      for(int iBin=1;iBin<=nBinsThisCat;iBin++) {
	catHists[iHist]->SetBinContent(iBin+offset,compHist->GetBinContent(iBin));
	
      }
    }
    offset+=nBinsThisCat;
  }

  // now vary only the corresponding bin
  offset=0;
  for(int iCat=0;iCat<catNames->size();iCat++) {
    //std::cout << iCat << "   " << offset << std::endl;
    TH2F* hist = (TH2F*)file->Get( "data_"+sms_pt+catNames->at(iCat)+"_"+postfix );
    std::unique_ptr<TH1F> compHist = compressHistogram(*hist,useVarBinning,1,1,&binningMap[catNames->at(iCat)]);
    //std::cout << histName << " " << postfix << " " << sms_pt << " " << iCat << "    " << dummy->Integral() << " --- " <<  compHist->Integral() << std::endl;
    int nBinsThisCat=nX*nY;
    if(useVarBinning) {      
      int nXBins = binningMap[ catNames->at(iCat) ].size();
      nBinsThisCat=nXBins*(nXBins-1)/2;
    }
    
    for(int iBin=1;iBin<=nBinsThisCat;iBin++) {
      catHists[iBin+offset-1]->SetBinContent(iBin+offset,compHist->GetBinContent(iBin));
    }
    offset+=nBinsThisCat;
  }
  //std::cout << ">>>>>bin 1:   " << catHist->GetBinContent(1) << std::endl;
  return catHists;
}

//a map between the bin of the input histogram and the bin of the output histogram
void CombinePrep::defineBinning(const TH2F& hist,std::vector<int>& xBinEdges, int minWidth,int minBin,float targetYield) {
  
  int maxBinX = hist.GetNbinsX();
  int maxBinY = hist.GetNbinsY();
  //std::cout << maxBinX << " " << maxBinY << std::endl;

  xBinEdges.push_back(maxBinX+1);

  int binLowX=maxBinX;
  int binLowY=1;
  
  TAxis *x = hist.GetXaxis();
  TAxis *y = hist.GetYaxis();

  while(binLowX>minBin) {
    while( hist.Integral(binLowX,maxBinX,binLowY,maxBinY) < targetYield  || maxBinX-binLowX<minWidth) {
      binLowX--;
      //std::cout << binLowX << "  " << hist.Integral(binLowX,maxBinX,binLowY,maxBinY) << "  " << maxBinX-binLowX << std::endl;
      if(binLowX==1 || binLowX==minBin) break;
    }
    //std::cout << binLowX << ":  " << hist.Integral(binLowX,maxBinX,binLowY,maxBinY) << std::endl;

    binLowY++;
    maxBinX=binLowX;
    xBinEdges.push_back(binLowX);    
  }

  std::cout << xBinEdges.size() << std::endl;
  for(int i=0;i<xBinEdges.size();i++) {
    std::cout << xBinEdges.at(i) << " ";
  }
  std::cout << std::endl;
  for(int i=0;i<xBinEdges.size()-1;i++) {
     for(int j=1;j<i+1;j++) {
       std::cout << x->GetBinLowEdge(xBinEdges.at(i+1)) << "--" <<x->GetBinLowEdge(xBinEdges.at(i)) << "   " << y->GetBinLowEdge(j) << "--"<<y->GetBinLowEdge(j+1) <<":    " << hist.Integral(xBinEdges.at(i+1),xBinEdges.at(i),j,j+1) << std::endl;
     }
     std::cout << x->GetBinLowEdge(xBinEdges.at(i+1)) << "--" <<x->GetBinLowEdge(xBinEdges.at(i)) << "   " << y->GetBinLowEdge(i+1) << "--"<<y->GetBinLowEdge(maxBinY+1) <<":    " << hist.Integral(xBinEdges.at(i+1),xBinEdges.at(i),i+1,maxBinY+1) << std::endl;
   }
}


void CombinePrep::defineExternalBinning(TString catName,std::vector<int>& xBinEdges) {
  externalBinning = true;

  if(binningMap.find(catName) == binningMap.end()) {
    binningMap.insert( std::pair<TString,std::vector<int>>(catName,{}) ); //setup the binning map
  }
  binningMap[catName].assign(xBinEdges.begin(),xBinEdges.end());
}


//takes in a single TH1F with correlated changes and returns a vector of single bin fluctuations.
//This allows easy conversion from a correlated to uncorrelated systematic
//histName should contain a %d where the index of the fluctuating bin will be placed 
std::vector<TH1F*> CombinePrep::uncorrelateHistogram(TH1F* fluc, TH1F* nom, TString histName) {

  std::vector<TH1F*> output;

  for(int i=0; i<= fluc->GetNbinsX(); i++) {
    TH1F* thisHist = (TH1F*)nom->Clone( Form(histName.Data(),i) );
    
    thisHist->SetBinContent(i, fluc->GetBinContent(i));
    output.push_back(thisHist);
  }
  return output;
}
