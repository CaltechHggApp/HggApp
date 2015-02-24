#include "makeSigRegionPredictionsFromTree.C"

#include "TString.h"

#include "SelectionStrings.h"

#include <iostream>
#include <cstdio>

using namespace SigRegionBinning;

void makeTTBarContributionTable () {

  TString mcFiles[] = {
    "~/raid4/TTj_RazorHgg.root",
    "~/raid4/MCBackground_DiPhotonJets_sherpa_plus_GJets_pythia_v4.root"
  };

  TFile ttFile(mcFiles[0]);
  TFile yyFile(mcFiles[1]);
  
  TTree* ttTree = (TTree*)ttFile.Get("SusyHggTree");
  TTree* yyTree = (TTree*)yyFile.Get("SusyHggTree");
  
  SelectionStrings sel;

  TFile ("TMP.root","RECREATE");
  TTree* ttSel = ttTree->CopyTree(sel.baseSelection);
  TTree* yySel = yyTree->CopyTree(sel.baseSelection);

  float target[5] = {358.7, 3.6, 6.1, 847.9, 1828.2};

  for(int i=0;i<5;i++) {
    std::vector<float> p_tt,p_yy;

    makeSigRegionPredictionsFromTree(ttSel->CopyTree(sel.mggSigRegion[i]+" && "+sel.boxDefs[i]),p_tt,static_cast<SigRegionBinning::BinningRegion>(i));
    makeSigRegionPredictionsFromTree(yySel->CopyTree(sel.mggSigRegion[i]+" && "+sel.boxDefs[i]),p_yy,static_cast<SigRegionBinning::BinningRegion>(i));

    std::vector<region> regs = getSigRegions(static_cast<SigRegionBinning::BinningRegion>(i));

    float sum = 0;
    for(int j=0;j<p_yy.size();j++) sum+=p_yy.at(j);

    std::cout << "SF: " << target[i]/sum << std::endl;

    for(int iReg=0; iReg<regs.size(); iReg++) {
      printf( " % 6.0f - % 6.0f & %0.2f - %0.2f & ",
	      regs.at(iReg).MR_min, regs.at(iReg).MR_max,
	      regs.at(iReg).Rsq_min, regs.at(iReg).Rsq_max);

      printf( " % 8.2f & % 8.2f \\\\\n",p_yy.at(iReg)*(target[i]-p_tt.at(iReg)*19.8/894)/sum,p_tt.at(iReg)*19.8/894);


    }
    std::cout << std::endl;
  }

}
