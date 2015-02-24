#ifndef SigRegionBinning_h
#define SigRegionBinning_h
#include <vector>
#include <algorithm>
#include <iostream>
#include "TRandom3.h"

#include "TH1D.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TMath.h"

namespace SigRegionBinning {
  enum BinningRegion {kHighPt=0, kHbb, kZbb, kHighRes, kLowRes};

  const int nBoxes=5;

  TString getRegionName(BinningRegion reg) {
    switch(reg) {
    case kHighPt:
      return "HighPt";
      
    case kHbb:
      return "Hbb";

    case kZbb:
      return "Zbb";

    case kHighRes:
      return "HighRes";

    case kLowRes:
      return "LowRes";
    }
    return "";
  }

  struct region {
    float MR_min;
    float MR_max;
    float Rsq_min;
    float Rsq_max;
    region(float minMr,float maxMr, float minRsq, float maxRsq) {
      MR_min=minMr; MR_max=maxMr; Rsq_min=minRsq; Rsq_max=maxRsq;
    }
  };
  
  //build the full list of regions from the MR edge list
  void buildRegionsFromEdges(std::vector<float>& MRedges, std::vector<region>& regions, float RsqStep=0.05) {
    std::sort(MRedges.begin(),MRedges.end()); //sort low to high
    for(int iMRbin=0;iMRbin<MRedges.size()-1;iMRbin++) {
      for(int iRsqBin=0;iRsqBin<MRedges.size()-1-iMRbin; iRsqBin++) {
	if(iRsqBin<MRedges.size()-iMRbin-2) { //not at the end in Rsq
	  regions.push_back(region(MRedges.at(iMRbin),MRedges.at(iMRbin+1),iRsqBin*RsqStep,(iRsqBin+1)*RsqStep));
	}
	else {
	  regions.push_back(region(MRedges.at(iMRbin),MRedges.at(iMRbin+1),iRsqBin*RsqStep,1.));	
	}
      }
    }
  }
  
  void printRegions(const std::vector<region>& regions) {
    std::vector<region>::const_iterator it = regions.begin();
    for(; it!=regions.end(); it++) {
      std::cout << it->MR_min << " - " << it->MR_max << "    " << it->Rsq_min << " - " << it->Rsq_max << std::endl;
    }
  }

  std::vector<float> getMREdges(BinningRegion r) {
    std::vector<float> edges;
    switch(r) {
    case kHighPt:
      edges.push_back(150);
      edges.push_back(200);
      edges.push_back(300);
      edges.push_back(500);
      edges.push_back(1600);
      edges.push_back(3000);
      break;
      
    case kHbb:
      edges.push_back(150);
      edges.push_back(300);
      edges.push_back(3000);
      
      break;
      
    case kZbb:
      edges.push_back(150);
      edges.push_back(450);
      edges.push_back(3000);
      break;
      
    case kHighRes:
      edges.push_back(150);
      edges.push_back(250);
      edges.push_back(400);
      edges.push_back(1400);
      edges.push_back(3000);
      break;
      
    case kLowRes:
      edges.push_back(150);
      edges.push_back(200);
      edges.push_back(250);
      edges.push_back(400);
      edges.push_back(1200);
      edges.push_back(3000);
      break;
    }
    return edges;
  }
  
  std::vector<region> getSigRegions(BinningRegion r) {
    std::vector<region> reg;
    std::vector<float> mrEdges = getMREdges(r);
    buildRegionsFromEdges( mrEdges , reg);
    
    //printRegions(reg);
    
    return reg;
  }


  TH2F* makeHistogram(BinningRegion r, TString name) {
    std::vector<float> vecMRedges = getMREdges(r);    
    float MRedges[10];
    float Rsqedges[10];

    std::copy(vecMRedges.begin(),vecMRedges.end(),MRedges);
    for(int i=0; i<=vecMRedges.size()-2; i++) {
      Rsqedges[i] = 0.05*i;
    }
    Rsqedges[vecMRedges.size()-1]=1.0;

    TH2F* hist = new TH2F(name,"",vecMRedges.size()-1,MRedges,vecMRedges.size()-1,Rsqedges);
    hist->SetXTitle("M_{R} (GeV)");
    hist->SetYTitle("R^{2}");

    return hist;
  }

  void formatSigRegionPlot(TH2F* hist) {
    const int nX = hist->GetNbinsX();
    const int nY = hist->GetNbinsY();

    for(int iXbin=1; iXbin<=nX; iXbin++) {
      float topContents = hist->GetBinContent(iXbin,nY-(iXbin-1));
      for(int iYbin=nY-(iXbin-1)+1; iYbin<=nY; iYbin++) {
	hist->SetBinContent(iXbin,iYbin,topContents);
      }
    }
  }

  TH1D* make1DProj(TH2F* hist,bool addError=false, float sf=1) {
    const int nX = hist->GetNbinsX();
    const int nY = hist->GetNbinsY();

    const int n1D = nX*(nX+1)/2;

    TH1D* h1D = new TH1D(TString(hist->GetName())+"_1D","",n1D,-0.5,n1D-0.5);
    TAxis* x = h1D->GetXaxis();
    h1D->SetYTitle("Events / Bin");

    int iBin=0;
    for(int iXbin=1; iXbin<=nX; iXbin++) {
      for(int iYbin=1; iYbin<nY-(iXbin-1)+1; iYbin++) {
	float content = hist->GetBinContent(iXbin,iYbin);
	if(iYbin==nY-(iXbin-1)) {
	  for(int iYIntBin = nY-(iXbin-1)+1; iYIntBin<=nY;iYIntBin++) {
	    content+=hist->GetBinContent(iXbin,iYIntBin);
	    std::cout << content << std::endl;
	  }
	}
	h1D->SetBinContent(++iBin,content);
	if(addError) {
	  h1D->SetBinError(iBin,TMath::Sqrt(content*(1+sf)));
	}
	x->SetBinLabel(iBin,Form("%0.2f #leq R^{2} < %0.2f  %0.0f #leq M_{R} < %0.0f",
				 hist->GetYaxis()->GetBinLowEdge(iYbin),
				 (iYbin==nY-(iXbin-1)?1.00: hist->GetYaxis()->GetBinLowEdge(iYbin+1)),
				 hist->GetXaxis()->GetBinLowEdge(iXbin),
				 hist->GetXaxis()->GetBinLowEdge(iXbin+1)
				 ));	

				 
      }
    }
    x->SetLabelSize(0.03);
    return h1D;
  }

  float getExpPer(TRandom3& rng, std::vector<float> expErr) {
    float expPer=0;
    for(std::vector<float>::const_iterator expIt=expErr.begin();
	expIt != expErr.end(); expIt++) {
      expPer+=rng.Gaus(0,*expIt);
    }
    return expPer;
  }

  //expErr should be a list of % errors for the expected background
  float pValue(float obs, float exp, std::vector<float> expErr) {
    //if(obs==0 && exp<0.50) return 1;
    TRandom3 rng;

    size_t count=0;
    float diff = fabs(obs-exp);

    size_t nToy=100;
    while(count<100) {
      nToy*=10;
      count=0;
      for(int i=0;i<nToy;i++) {
	float thisExpMean = exp*getExpPer(rng,expErr);
	float thisExp = rng.Poisson(thisExpMean);
	if( fabs(thisExp-obs) >= diff ) count++;
      }
    }
    return float(count)/nToy;
  }

  float pValue(float obs, float exp, float expErr) {
    std::vector<float> vecExpErr(1,expErr);
    return pValue(obs,exp,vecExpErr);
  }

};

#endif
