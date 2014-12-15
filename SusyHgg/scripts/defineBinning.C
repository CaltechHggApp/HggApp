#include "TH2F.h"
#include "TAxis.h"
#include <iostream>

//a map between the bin of the input histogram and the bin of the output histogram
void defineBinning(const TH2F& hist) { //,std::vector<float>& binEdges,float targetYield=10.) {
  
  int min_width=5;
  int min = 20;

  int maxBinX = hist.GetNbinsX();
  int maxBinY = hist.GetNbinsX();
  std::vector<int> binEdges;
  float targetYield=4.;
  binEdges.push_back(maxBinX+1);

  int binLowX=maxBinX;
  int binLowY=1;
  
  TAxis *x = hist.GetXaxis();
  TAxis *y = hist.GetYaxis();

  while(binLowX>min) {
    while( hist.Integral(binLowX,maxBinX,binLowY,maxBinY) < targetYield  || maxBinX-binLowX<min_width) {
      binLowX--;
      if(binLowX==1) break;
    }
    std::cout << binLowX << ":  " << hist.Integral(binLowX,maxBinX,binLowY,maxBinY) << std::endl;

    binLowY++;
    maxBinX=binLowX;
    binEdges.push_back(binLowX);    
  }


  for(int i=0;i<binEdges.size()-1;i++) {
    for(int j=1;j<i+1;j++) {
      std::cout << x->GetBinLowEdge(binEdges.at(i+1)) << "--" <<x->GetBinLowEdge(binEdges.at(i)) << "   " << y->GetBinLowEdge(j) << "--"<<y->GetBinLowEdge(j+1) <<":    " << hist.Integral(binEdges.at(i+1),binEdges.at(i),j,j+1) << std::endl;
    }
    std::cout << x->GetBinLowEdge(binEdges.at(i+1)) << "--" <<x->GetBinLowEdge(binEdges.at(i)) << "   " << y->GetBinLowEdge(i+1) << "--"<<y->GetBinLowEdge(maxBinY+1) <<":    " << hist.Integral(binEdges.at(i+1),binEdges.at(i),i+1,maxBinY+1) << std::endl;
  }
  

}
