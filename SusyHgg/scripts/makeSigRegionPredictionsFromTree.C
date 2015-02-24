#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <cstdio>
#include "SigRegionBinning.h"

using namespace SigRegionBinning;

void makeSigRegionPredictionsFromTree(TTree *tree, std::vector<float >& predictions, BinningRegion reg) {
  std::vector<region> regs = getSigRegions(reg);
  
  for(int i=0;i<regs.size();i++) {
    predictions.push_back( tree->GetEntries( Form( "MR>=%f && MR<%f && t1Rsq>=%f && t1Rsq<%f",
						  regs.at(i).MR_min,
						  regs.at(i).MR_max,
						  regs.at(i).Rsq_min,
						  regs.at(i).Rsq_max )
					    )
			  );

  }

}

void printSigRegionPredictionsFromTree(TTree *tree, BinningRegion reg,float norm) {
  std::vector<float > predictions;
  std::vector<region> regs = getSigRegions(reg);
  makeSigRegionPredictionsFromTree(tree,predictions, reg);
  for(int i=0;i<regs.size();i++) {
    printf( "% 7.0f - % 7.0f & %0.2f - %0.2f & % 10.3f \\\\\n",
	    regs.at(i).MR_min, regs.at(i).MR_max,
	    regs.at(i).Rsq_min, regs.at(i).Rsq_max,
	    predictions.at(i)*norm);
  }
}

void printSigRegionPredictionsFromManyTrees(TTree **trees,int nTrees, BinningRegion reg,float norm) {
  std::vector<region> regs = getSigRegions(reg);
  std::vector< std::vector<float> > predictions(nTrees);
  for(int iTree=0; iTree<nTrees; iTree++) {
    makeSigRegionPredictionsFromTree(trees[iTree],predictions.at(iTree), reg);
  }

  for(int i=0;i<regs.size();i++) {
    printf( "% 7.0f - % 7.0f & %0.2f - %0.2f ",
	    regs.at(i).MR_min, regs.at(i).MR_max,
	    regs.at(i).Rsq_min, regs.at(i).Rsq_max);

    for(int iTree=0;iTree<nTrees;iTree++) {
      printf("& %10.3f ",predictions.at(iTree).at(i));
    }
    printf("\\\\\n");
  }
}
