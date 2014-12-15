#include <vector>
#include <algorithm>
#include <iostream>

enum BinningRegion {kHighPt, kHbb, kZbb, kHighRes, kLowRes};

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

std::vector<region> getSigRegions(BinningRegion r) {
  std::vector<region> reg;

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

  buildRegionsFromEdges( edges, reg);

  //printRegions(reg);

  return reg;
}
