#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TList.h"
#include "TRandom3.h"
#include "TMath.h"
#include "RooStats/RooStatsUtils.h"

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <map>

#include "SigRegionBinning.h"
#include "assert.h"
#include <stdio.h>
#include "profileNoHistMinuit.C"

using namespace std;
using namespace SigRegionBinning;

#define NEWPROC 0

#define MAX(a,b) ((a)>(b) ? (a) : (b))

typedef vector<vector<pair<float,float> > > systMap;

void convertHistVector(const TH1F* hist, vector<float>& vec) {
  for(int i=1;i<=hist->GetNbinsX();i++) {vec.push_back(hist->GetBinContent(i));}
}

void convertHistPairVector(const TH1F* hist1, const TH1F* hist2, vector<pair<float,float> >& vec) {
  for(int i=1;i<=hist1->GetNbinsX();i++) {vec.push_back(make_pair(hist1->GetBinContent(i),hist2->GetBinContent(i)));}
}

void pruneSystematics(float nom, vector<pair<float,float> >& systs, float threshold=0.1);

pair<float,float> addSystsQuad(float nom, vector<pair<float,float> >& systs);

float throwToyMean(TRandom3& rng, float nom, vector<pair<float,float> >& systs);

float getSigMinuit(float obs, std::pair<float,float> sideband, std::pair<float,float> SFHigh, std::pair<float,float> SFLow,
	   std::pair<float,float> higgsBkg);

void makeUnblindTable_vProfile(TString dir="./", bool BLIND=true,bool fullTex=false) {
  TFile dataFile(dir+"/data/data.root");
  if(!dataFile.IsOpen()) {
    cout << "Cannot open " << dir+"/data.root" << endl;
    return;
  }

  TH1F* obs;
  if( BLIND ) {
    obs = (TH1F*)dataFile.Get("bkg");
    cout << "***************\n"
	 << "*   BLINDED   *\n"
	 << "***************\n";

  } else {
    obs = (TH1F*)dataFile.Get("data_obs");
  }

  TH1F* scaleFactorHist = (TH1F*)dataFile.Get("scaleFactors");
  if(scaleFactorHist==0) {
    std::cout << "OLD DATA FORMAT DETECTED!!!\nABORTING" << std::endl;
    return;
  }


    if(obs==0) {
      std::cout << "cannot find histogram data_obs" << std::endl;
      return;
    }

    TH1F* bkgStatHist     = (TH1F*)dataFile.Get("bkg_statistics");
    TH1F* bkgStatHighHist = (TH1F*)dataFile.Get("bkg_bkgShapeHigh");
    TH1F* bkgStatLowHist = (TH1F*)dataFile.Get("bkg_bkgShapeLow");
    if(bkgStatHist==0 || bkgStatHighHist==0 || bkgStatLowHist==0) {
      std::cout << "Cannot get background statistics" << std::endl;
      return;
    }

    float totalStat = bkgStatHist->Integral();
    float totalHigh = bkgStatHighHist->Integral();
    float totalLow = bkgStatLowHist->Integral();

    const int Nbkg=5;
    TString bkgNames[Nbkg] = {"bkg","ggH", "vbfH", "wzH", "ttH"};
    //const int Nbkg=1;
    //TString bkgNames[Nbkg] = {"bkg"};

  TList * keys = dataFile.GetListOfKeys();

  //bin
  vector<float> obsVec;
  convertHistVector(obs,obsVec);

  vector<float>   scaleFactors;
  vector<float>   scaleFactorsError;
  convertHistVector(scaleFactorHist,scaleFactors);

  //pair(SF,SFE)
  std::vector<std::pair<float,float>> scaleFactorHigh;
  std::vector<std::pair<float,float>> scaleFactorLow;

  vector<float>   bkgStatistics;
  convertHistVector(bkgStatHist,bkgStatistics);

  vector<std::pair<float,float>>   bkgStatisticsHighLow;
  for(int i=1;i<=bkgStatHighHist->GetNbinsX();i++) {
    bkgStatisticsHighLow.push_back(std::make_pair(bkgStatHighHist->GetBinContent(i),bkgStatLowHist->GetBinContent(i)));
  }

  //wow, this is ugly...
  vector<float> totalStatistics; // the total statistics of the box this bin belongs to
  vector<std::pair<float,float>> totalStatisticsHighLow; // the total statistics of the box this bin belongs to
  int index=0;
  for(int iReg=0; iReg<nBoxes; iReg++) {
    vector<region> region = getSigRegions((BinningRegion)iReg);

    std::pair<float,float> tot = make_pair(0,0);
    for(int i=0;i<region.size();i++) {
      int iC = ConvertToCombineBin((BinningRegion)iReg,i);
      tot.first+=bkgStatisticsHighLow.at(i).first;
      tot.second+=bkgStatisticsHighLow.at(i).second;
    }
    cout << tot.first << endl;
    for(int i=0;i<region.size();i++) {
      totalStatisticsHighLow.push_back(tot);
      totalStatistics.push_back(tot.first+tot.second);
    }

  }

  //bin<bkg>
  vector<vector<float> >   bkgNominal(obsVec.size());
  //bin<bkg<systematic>>
  vector<systMap>          bkgSyst(obsVec.size(),systMap(Nbkg));
  vector<vector<vector<TString> > > bkgSystNames(obsVec.size(),vector<vector<TString> >(Nbkg));

  

  //per signal region bin
  for(int iBin=0; iBin<obsVec.size(); iBin++) { 
    float shift = (NEWPROC?scaleFactors.at(iBin):0);
    
    TH1F* SFErrHist = (TH1F*)dataFile.Get("bkg_fitUp");
    
    //per background component
    for(int iBkg=0; iBkg<Nbkg; iBkg++) {
      TH1F* nomH = (TH1F*)dataFile.Get(bkgNames[iBkg]);
      if(nomH==0) {
	std::cout << "cannot find histogram " << bkgNames[iBkg] << std::endl;
	return;
      }
      bkgNominal.at(iBin).push_back(nomH->GetBinContent(iBin+1)+shift);

      if(iBkg==0) {
	scaleFactorsError.push_back( scaleFactors.at(iBin)* (SFErrHist->GetBinContent(iBin+1)/bkgNominal.at(iBin).at(0)-1) );
	//std::cout << scaleFactors.at(iBin) << " +- " << scaleFactorsError.at(iBin) << std::endl;

	
	float sfH = scaleFactors.at(iBin)/totalStatisticsHighLow.at(iBin).first*totalStatistics.at(iBin);
	float sfL = scaleFactors.at(iBin)/totalStatisticsHighLow.at(iBin).second*totalStatistics.at(iBin);
	
	float sfHE = sfH/scaleFactors.at(iBin)*TMath::Sqrt( TMath::Power(scaleFactorsError.at(iBin),2)
				  - 1/(totalStatistics.at(iBin))
				  + 1/totalStatisticsHighLow.at(iBin).first);
	float sfLE = sfL/scaleFactors.at(iBin)*TMath::Sqrt( TMath::Power(scaleFactorsError.at(iBin),2)
				  - 1/(totalStatistics.at(iBin))
				  + 1/totalStatisticsHighLow.at(iBin).second);
	
	//std::cout << sfH << " +- " << sfHE << std::endl;
	//std::cout << sfL << " +- " << sfLE << std::endl;
	//std::cout << std::endl;
	
	scaleFactorHigh.push_back(std::make_pair(sfH,sfHE));
	scaleFactorLow.push_back(std::make_pair(sfL,sfLE));
      }

      //bkgSyst.at(iBin).reserve(Nbkg);
      

      //per systematic
      for(int iSyst=0; iSyst<keys->GetEntries(); iSyst++) {
	const char* name = keys->At(iSyst)->GetName();
	//no match or its just the nominal
	if(strncmp(bkgNames[iBkg].Data(),name,bkgNames[iBkg].Length())!=0 
	   || bkgNames[iBkg].Length()==strlen(name)) {
	  //delete name;
	  continue;
	}
	
	TString tsName(name);
	//delete name;      

	//only run for the Up systematic, we'll process both together to form the pair
	if(! tsName.EndsWith("Up")) {
	  continue;
	}
	
	systMap systematics;
	
	TH1F* up   = (TH1F*)dataFile.Get(tsName);
	TH1F* down = (TH1F*)dataFile.Get(tsName(0,tsName.Length()-2)+"Down");
	if(up==0) {
	  cout << "could not find histogram named " << tsName << endl;
	  return;
	}
	if(down==0) {
	  cout << "could not find histogram named " << tsName(0,tsName.Length()-2)+"Down" << endl;
	  return;
	}
	if(up->GetBinContent(iBin+1)==down->GetBinContent(iBin+1)) continue;
	//deal with the sideband statistcs systematic differently
	if(tsName.Contains("sidebandStatistics") || tsName.Contains("bkgShape") || tsName.Contains("bkg_fitUp")) {
	  continue;
	  float N;
	  if( up->GetBinContent(iBin+1)/bkgNominal.at(iBin).at(iBkg)==1 ) N=1;
	  else N = TMath::Power(1./(up->GetBinContent(iBin+1)/(bkgNominal.at(iBin).at(iBkg)-shift)-1),2); //compute the statistics out of the card
	  float sf = (up->GetBinContent(iBin+1)+shift-bkgNominal.at(iBin).at(iBkg))/sqrt(N);
	  bkgStatistics.push_back( N );
	  //scaleFactors.push_back( sf);
	  std::cout << obsVec.at(iBin) << "  " << bkgNominal.at(iBin).at(iBkg) << "   " << up->GetBinContent(iBin+1) << "  " << N << "  " << sf << std::endl;
	  bkgSyst.at(iBin).at(iBkg).push_back( make_pair(bkgNominal.at(iBin).at(iBkg)+sqrt(N+0.76)*sf+shift,bkgNominal.at(iBin).at(iBkg)-sqrt(N+0.76)*sf+shift) );
	  bkgSystNames.at(iBin).at(iBkg).push_back(tsName);
	}else{
	  bkgSyst.at(iBin).at(iBkg).push_back( make_pair(up->GetBinContent(iBin+1)+shift,down->GetBinContent(iBin+1)+shift));
	  bkgSystNames.at(iBin).at(iBkg).push_back(tsName);
	}
      }
    }
  }

  for(int iReg=0; iReg<nBoxes; iReg++) {
    vector<region> region = getSigRegions((BinningRegion)iReg);
    cout << endl <<endl << getRegionName( (BinningRegion)iReg ) << endl << endl;

    if(fullTex) {
      std::cout << "\\begin{table}\n\\begin{tabular}{|cc|c|c|cc|}" << std::endl
		<< "\\hline\n"
		<< "$M_R$ region & $R^2$ region & observed events & expected background & p-value & significance ($\\sigma$) \\\\\n" 
		<< "\\hline" << std::endl;
	}
    
    for(int iC=0;iC<region.size();iC++) {
      int i = ConvertToCombineBin((BinningRegion)iReg,iC);
      float bkgTot=0;
      for(int iBkg=0; iBkg < Nbkg; iBkg++) bkgTot+=bkgNominal.at(i).at(iBkg);
      
      if(obsVec.at(i)==5 && false) {
	for(int iBkg=0; iBkg<Nbkg; iBkg++) {
	  std::cout << bkgNames[iBkg] << ":  " << bkgNominal.at(i).at(iBkg) << std::endl;
	  for(int iSyst=0; iSyst< bkgSyst.at(i).at(iBkg).size(); iSyst++) {
	    std::cout << "\t" << bkgSystNames.at(i).at(iBkg).at(iSyst) << "\t:" << bkgSyst.at(i).at(iBkg).at(iSyst).second << " -- " << bkgSyst.at(i).at(iBkg).at(iSyst).first << std::endl;
	  }
	}
      }
      
      float higgsTot=0;
      pair<float,float> err = make_pair(0,0);
      for(int iBkg=1; iBkg < Nbkg; iBkg++) { 
	//std::cout << bkgNominal.at(i).at(iBkg);
	higgsTot+=bkgNominal.at(i).at(iBkg);
	pair<float,float> thisErr = addSystsQuad(bkgNominal.at(i).at(iBkg),bkgSyst.at(i).at(iBkg));
	err.first+=thisErr.first;
	err.second+=thisErr.second;
      }
      //err.first+=(bkgStatistics.at(i)+1)*scaleFactors.at(i)*scaleFactors.at(i);
      //err.second+=(bkgStatistics.at(i)+1)*scaleFactors.at(i)*scaleFactors.at(i);

      pair<float,float> dispErr = err; // the error to display (add in the statistical and sideband errors)

      dispErr.first = dispErr.first + scaleFactors.at(i)*scaleFactors.at(i)*bkgStatistics.at(i); //add in the background stat error (in quadrature)
      dispErr.second = dispErr.second + scaleFactors.at(i)*scaleFactors.at(i)*bkgStatistics.at(i); //add in the background stat error (in quadrature)

      std::pair<float,float> sideband = bkgStatisticsHighLow.at(i);

      if( fabs(sideband.first-sideband.second) < sqrt(bkgStatistics.at(i)) ) {
	//sideband error small, use statistical error
	dispErr.first = dispErr.first + 0.5*scaleFactors.at(i)*scaleFactors.at(i)*bkgStatistics.at(i); //add in the background stat error (in quadrature)
	dispErr.second = dispErr.second + 0.5*scaleFactors.at(i)*scaleFactors.at(i)*bkgStatistics.at(i); //add in the background stat error (in quadrature)
      } else {
	float diff = fabs(sideband.first-sideband.second);
	dispErr.first = dispErr.first + 0.5*scaleFactors.at(i)*scaleFactors.at(i)*diff; //add in the background shape systematics
	dispErr.second = dispErr.second + 0.5*scaleFactors.at(i)*scaleFactors.at(i)*diff;
      }

      // print \pm if the up/down errors are the same
      if( fabs(sqrt(dispErr.first) - sqrt(dispErr.second)) < 0.01 ) {
	printf( "% 6.0f - % 6.0f & %0.2f - %0.2f & % 4.0f & $% 4.1f \\pm %0.2f$ ",
		region.at(iC).MR_min, region.at(iC).MR_max,
	region.at(iC).Rsq_min, region.at(iC).Rsq_max,	    
		obsVec.at(i), bkgTot,sqrt(dispErr.first)); 
      }else{
	printf( "% 6.0f - % 6.0f & %0.2f - %0.2f & % 4.0f & $% 4.1f^{+%0.2f}_{-%0.2f}$ ",
		region.at(iC).MR_min, region.at(iC).MR_max,
		region.at(iC).Rsq_min, region.at(iC).Rsq_max,	    
		obsVec.at(i), bkgTot,sqrt(dispErr.first),sqrt(dispErr.second)); 
      }
      cout << flush;

      float obs = obsVec.at(i);
      std::pair<float,float> SFHigh   = scaleFactorHigh.at(i);
      std::pair<float,float> SFLow   = scaleFactorLow.at(i);
      std::pair<float,float> higgs = make_pair(higgsTot,sqrt(err.first));

      

      float sig = getSigMinuit(obs,sideband,SFHigh,SFLow,higgs);
      if(isinf(sig)) return;

      float pv = RooStats::SignificanceToPValue(sig);//pval(obsVec.at(i),bkgNominal.at(i),bkgSyst.at(i));
      //float sig = fabs(TMath::NormQuantile(pv/2));
      if(obsVec.at(i) < bkgTot && fabs(sig)>0.005) sig*=-1;
      printf("& %0.3f & %0.1f \\\\\n",pv, sig);
    }
    if(fullTex) {
      std::cout << "\\hline\n\\end{tabular}" << std::endl
		<< "\\caption{Number of Events observed in the signal region compared to expected background in the {\\bf " << getRegionName( (BinningRegion)iReg ) << "} box.}\n"
		<< "\\label{tab:ObsPV_" << getRegionName( (BinningRegion)iReg ) << "}\n"
		<< "\\end{table}" << std::endl;
    }
    //break;
  }
  
  
}


pair<float,float> addSystsQuad(float nom, vector<pair<float,float> >& systs) {
  pair<float,float> sumQuad = make_pair(0.0,0.0);

  for(int i=0; i<systs.size(); i++) {
    sumQuad.first  += TMath::Power(nom - systs.at(i).first,2);
    sumQuad.second += TMath::Power(nom - systs.at(i).second,2);
  }
  return sumQuad;
}


float throwToyMean(TRandom3& rng, float nom, vector<pair<float,float> >& systs) {
  float mean = nom;
  //#if 0
  for(int iSyst=0; iSyst<systs.size(); iSyst++) {
    bool up = (rng.Uniform(1) < 0.5); // whether this toy fluctuates up or down
    float thisSyst = fabs(nom - (up ? systs.at(iSyst).first : systs.at(iSyst).second));
    float delta = fabs(rng.Gaus(0,thisSyst));
    if(up) mean+=delta;
    else   mean-=delta;

    //std::cout << "\t\t" << up << "  " << thisSyst << "  " << delta << "  " << mean << endl;
  }
  //#endif
  return mean;
}


float pval(float obs, vector<float> mean_exp, systMap& systs) {
  TRandom3 rng(0);
  size_t count = 0;
  float totalBkg = 0;
  for(int i=0;i<mean_exp.size();i++) totalBkg+=mean_exp.at(i);
  float diff = fabs(obs-totalBkg);

  size_t nToy=100;
  while(count<1000) {
    nToy*=10;
    count=0;
    for(size_t i=0; i<nToy; i++) {
      float toyMean=0;
      for(int iBkg=0;iBkg<mean_exp.size();iBkg++) {
	toyMean += throwToyMean(rng,mean_exp.at(iBkg),systs.at(iBkg));
	//if(i%100==0) { cout << mean_exp.at(iBkg) << "  " << toyMean <<endl; }
      }
      if( fabs(rng.Poisson(MAX(toyMean,0))-totalBkg)>=diff) count++;
    }
  }
  return float(count)/nToy;
}


void pruneSystematics(float nom, vector<pair<float,float> >& systs, float threshold) {
  float largestSyst=0;
  vector<pair<float,float> >::iterator systs_it=systs.begin();
  for(; systs_it!=systs.end(); systs_it++) {
    largestSyst = MAX(largestSyst, fabs(systs_it->first-nom));
    largestSyst = MAX(largestSyst, fabs(systs_it->second-nom));
  }

  systs_it=systs.begin();
  while(systs_it!=systs.end()) {
    if(fabs(systs_it->first-nom)<threshold*largestSyst && fabs(systs_it->second-nom)<threshold*largestSyst) {
      systs_it = systs.erase(systs_it);
    } else {
      systs_it++;
    }    
  }
}

float getSigMinuit(float obs, std::pair<float,float> sideband, std::pair<float,float> SFHigh, std::pair<float,float> SFLow,
		   std::pair<float,float> higgsBkg) {


  params.obs = obs;
  params.Nupper = sideband.second;
  params.Nlower = sideband.first;
  //params.Nupper = sideband.first;
  //params.Nlower = sideband.second;
  params.UpperSF = SFHigh.first;
  params.UpperSFe = SFHigh.second;
  params.LowerSF = SFLow.first;
  params.LowerSFe = SFLow.second;
  params.Higgs = higgsBkg.first;
  params.HiggsErr = higgsBkg.second;

  std::pair<float,float> prof = profileNoHistMinuit();

  return sqrt(prof.first);
}
