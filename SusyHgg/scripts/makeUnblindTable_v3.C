#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TList.h"
#include "TRandom3.h"
#include "TMath.h"

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <map>

#include "SigRegionBinning.h"
#include "assert.h"

using namespace std;
using namespace SigRegionBinning;



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

float pval(float obs, vector<float> mean_exp, systMap& systs);

void makeUnblindTable_v2(TString dir="./", bool BLIND=true,bool fullTex=false) {
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



    if(obs==0) {
      std::cout << "cannot find histogram data_obs" << std::endl;
      return;
    }

  const int Nbkg=5;
  TString bkgNames[Nbkg] = {"bkg","ggH", "vbfH", "wzH", "ttH"};

  TList * keys = dataFile.GetListOfKeys();

  //bin
  vector<float> obsVec;
  convertHistVector(obs,obsVec);
  //bin<bkg>
  vector<vector<float> >   bkgNominal(obsVec.size());
  //bin<bkg<systematic>>
  vector<systMap>          bkgSyst(obsVec.size(),systMap(Nbkg));
  vector<vector<vector<TString> > > bkgSystNames(obsVec.size(),vector<vector<TString> >(Nbkg));
  //per signal region bin
  for(int iBin=0; iBin<obsVec.size(); iBin++) { 
    
    
    //per background component
    for(int iBkg=0; iBkg<Nbkg; iBkg++) {
      TH1F* nomH = (TH1F*)dataFile.Get(bkgNames[iBkg]);
      if(nomH==0) {
	std::cout << "cannot find histogram " << bkgNames[iBkg] << std::endl;
	return;
      }
      bkgNominal.at(iBin).push_back(nomH->GetBinContent(iBin+1));

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
	bkgSyst.at(iBin).at(iBkg).push_back( make_pair(up->GetBinContent(iBin+1),down->GetBinContent(iBin+1)));
	bkgSystNames.at(iBin).at(iBkg).push_back(tsName);
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
      
      pair<float,float> err = make_pair(0,0);
      for(int iBkg=0; iBkg < Nbkg; iBkg++) { 
	pair<float,float> thisErr = addSystsQuad(bkgNominal.at(i).at(iBkg),bkgSyst.at(i).at(iBkg));
	err.first+=thisErr.first;
	err.second+=thisErr.second;
      }

      printf( "% 6.0f - % 6.0f & %0.2f - %0.2f & % 4.0f & $% 4.2f^{+%0.3f}_{-%0.3f}$ ",
	      region.at(iC).MR_min, region.at(iC).MR_max,
	      region.at(iC).Rsq_min, region.at(iC).Rsq_max,	    
	      obsVec.at(i), bkgTot,sqrt(err.first),sqrt(err.second)); 
      cout << flush;
      float pv = pval(obsVec.at(i),bkgNominal.at(i),bkgSyst.at(i));
      printf("& %0.6f & %0.2f \\\\\n",pv, -1*TMath::NormQuantile(pv/2));
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
    sumQuad.first  += TMath::Power(nom-systs.at(i).first,2);
    sumQuad.second += TMath::Power(nom-systs.at(i).second,2);
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
  while(count<100) {
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
