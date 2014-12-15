#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>

#include "SigRegionBinning.h"

void makeSigRegionPredictionsFromTree(TTree *tree, std::vector<float >& predictions, BinningRegion reg, bool print=false) {
  std::vector<region> regions = getSigRegions(reg);

  for(int i=0;i<regions.size();i++) {
    predictions.push_back( tree->GetEntries( Form( "MR>=%f && MR<%f && t1Rsq>=%f && t1Rsq<%f",
						  regions.at(i).MR_min,
						  regions.at(i).MR_max,
						  regions.at(i).Rsq_min,
						  regions.at(i).Rsq_max )
					    )
			  );
    if(print) {
      std::cout << regions.at(i).MR_min << " - " << regions.at(i).MR_max << "    "
		<< regions.at(i).Rsq_min << " - " << regions.at(i).Rsq_max << "    "
		<< predictions.at(i) << std::endl;
    }

  }

}
