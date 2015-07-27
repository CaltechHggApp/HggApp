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
#include "profileSimpleNoHistMinuit.C"

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

float getSigMinuit( float obs, float sideband, std::pair<float,float> SF,
		    std::pair<float,float> higgsBkg, float bkgShapeErr );

TH1D* getDeltaLogLikelihood(float obs, float sideband, std::pair<float,float> SF,
			    std::pair<float,float> higgsBkg,float bkgShapeErr, bool _profileNobs = false );

pair<float, float> findOneSigma( TH1D* _nll );

void makeUnblindTable_vProfile( TString dir="./", bool BLIND=true, bool fullTex=false ) 
{
  //-------------------------------------
  //open runSusyHggDir/data/dat.root
  //contains all bins including all categories
  //-------------------------------------
  TFile dataFile(dir+"/data/data.root");
  if(!dataFile.IsOpen()) 
    {
      cout << "Cannot open " << dir+"/data.root" << endl;
      return;
    }

  //-----------------------------------------------------
  //getting data_obs histo, containing all 46 MR-Rsq bins
  //-----------------------------------------------------
  TH1F* obs;
  if( BLIND ) {
    obs = (TH1F*)dataFile.Get("bkg");
    cout << "***************\n"
	 << "*   BLINDED   *\n"
	 << "***************\n";

  } else {
    obs = (TH1F*)dataFile.Get("data_obs");
  }

  //-------------------------------------------------------
  //getting scale factor histogram
  //filled with boxes scale factors 
  //(
  // there are actually 46 bin in this histo.
  // 15-HighPt, 3-Hbb, 3-Zbb, 10-HighRes, 15-LowRes,
  // each 15-HighPt scale factor are numerical very close,
  // but not identical, why is that?
  //)
  //-------------------------------------------------------
  TH1F* scaleFactorHist = (TH1F*)dataFile.Get("scaleFactors");
  if( scaleFactorHist==0 ) {
    std::cout << "OLD DATA FORMAT DETECTED!!!\nABORTING" << std::endl;
    return;
  }
  
  
  if(obs==0) {
    std::cout << "cannot find histogram data_obs" << std::endl;
    return;
  }
  
  //------------------------------------------------------
  //getting sideband statistics histograms (46 MR-Rsq bins)
  //------------------------------------------------------
  TH1F* bkgStatHist     = (TH1F*)dataFile.Get("bkg_statistics");//this one should be the nsidebands events
  TH1F* bkgStatHighHist = (TH1F*)dataFile.Get("bkg_bkgShapeHigh");//this should be nhigh_sideband events 
  TH1F* bkgStatLowHist  = (TH1F*)dataFile.Get("bkg_bkgShapeLow");//this should be nlow_sideband events
  if( bkgStatHist == 0 || bkgStatHighHist == 0 || bkgStatLowHist == 0 ) 
    {
      std::cout << "Cannot get background statistics" << std::endl;
      return;
    }
  
  float totalStat = bkgStatHist->Integral();
  float totalHigh = bkgStatHighHist->Integral();
  float totalLow  = bkgStatLowHist->Integral();
  
  //-----------------------------------------------------
  //Defining combinatorial and SM-Higgs background labels
  //-----------------------------------------------------
  const int Nbkg=5;
  TString bkgNames[Nbkg] = {"bkg","ggH", "vbfH", "wzH", "ttH"};
    
  TList * keys = dataFile.GetListOfKeys();

  //-----------------------------------------------------------
  //converting histo (46 MR-Rsq bins) into a vector for obs
  //-----------------------------------------------------------
  vector<float> obsVec;
  convertHistVector(obs,obsVec);

  //--------------------------------------------------------------------
  //converting histo (46 MR-Rsq bins) into a vecor of SFs
  //--------------------------------------------------------------------
  vector<float>   scaleFactors;
  vector<float>   scaleFactorsError;
  convertHistVector(scaleFactorHist,scaleFactors);

  //------------------------------------------------------------------------------
  //define vector of pairs<SF,SFerror> for the high and low sideband scale factors
  //------------------------------------------------------------------------------
  std::vector<std::pair<float,float>> scaleFactorHigh;
  std::vector<std::pair<float,float>> scaleFactorLow;

  //-------------------------------------------------------------
  //conveting bkg statistics histo into a vector (46 MR-Rsq bins)
  //-------------------------------------------------------------
  vector<float>   bkgStatistics;
  convertHistVector(bkgStatHist,bkgStatistics);

  //--------------------------------------------------------------------------------------
  //creating vector of pairs<Nsideband_high, Nsideband_low> for each of the 46 MR-Rsq bins
  //--------------------------------------------------------------------------------------
  vector<std::pair<float,float>>   bkgStatisticsHighLow;
  for(int i=1;i<=bkgStatHighHist->GetNbinsX();i++) {
    bkgStatisticsHighLow.push_back(std::make_pair(bkgStatHighHist->GetBinContent(i),bkgStatLowHist->GetBinContent(i)));
  }

  for ( int k = 0; k < bkgStatistics.size(); k++ )
    {
      std::cout << "[DEBUG]: Nsideband: " << bkgStatistics.at(k) 
		<< " Nsideband-high: " << bkgStatisticsHighLow.at(k).first
		<< " Nsideband-low: " << bkgStatisticsHighLow.at(k).second
		<< std::endl;
    }

  //wow, this is ugly...
  vector<float> totalStatistics; // the total statistics of the box this bin belongs to
  vector<std::pair<float,float>> totalStatisticsHighLow; // the total statistics of the box this bin belongs to
  int index=0;
  //nBoxes,region define in namespace SigRegionBinning (see SigRegionBinning.h)
  for( int iReg = 0; iReg < nBoxes; iReg++ ) 
    {
      std::vector<region> region = getSigRegions((BinningRegion)iReg);
      std::pair<float,float> tot = make_pair(0,0);
      std::cout << "[INFO]: BoxName: " << getRegionName( (BinningRegion)iReg ) 
		<< " ,number of bins (signal regions): " << region.size() 
		<< std::endl;
      for( int i = 0; i < region.size(); i++ ) 
	{
	  int iC = ConvertToCombineBin((BinningRegion)iReg,i);
	  std::cout << "[DEBUG]: iC->" << iC << " bkgStatisticsHighLow: " << bkgStatisticsHighLow.at(iC).first
		    << ", " << bkgStatisticsHighLow.at(iC).second << std::endl;
	  tot.first  += bkgStatisticsHighLow.at(iC).first;
	  tot.second += bkgStatisticsHighLow.at(iC).second;
	}
      std::cout << "nsideband high: " << tot.first 
		<< " , nsideband low: " << tot.second 
		<< std::endl;
      for( int i = 0; i < region.size(); i++ ) 
	{
	  totalStatisticsHighLow.push_back(tot);
	  totalStatistics.push_back(tot.first+tot.second);
	}
    }
  
  for ( int k = 0; k < totalStatisticsHighLow.size(); k++ )
    {
      std::cout << "[DEBUG]: i->" << k << " HighLow Statistics: (" 
		<< totalStatisticsHighLow.at(k).first << "," << totalStatisticsHighLow.at(k).second 
		<< "); Full Statistics: " << totalStatistics.at(k) << std::endl;
    }

  //bin<bkg>
  vector<vector<float> >   bkgNominal(obsVec.size());
  std::cout << "[DEBUG]: bkgNominal size: " << bkgNominal.size() 
	    << " subvector size: " << bkgNominal.at(0).size() << std::endl;
  //bin<bkg<systematic>>
  vector<systMap>          bkgSyst(obsVec.size(),systMap(Nbkg));
  vector<vector<vector<TString> > > bkgSystNames(obsVec.size(),vector<vector<TString> >(Nbkg));

  
  //-----------------------------------------
  //looping over 46 MR-Rsq bins
  //-----------------------------------------
  for( int iBin = 0; iBin < obsVec.size(); iBin++ ) 
    { 
      float shift = (NEWPROC?scaleFactors.at(iBin):0);//in this case shift is always zero if NEWPROC is zero
      TH1F* SFErrHist = (TH1F*)dataFile.Get("bkg_fitUp");
      
      //--------------------------
      //loop over bkg types (non-resonant, ggH, vbfH, wzH, ttH)
      //--------------------------
      for( int iBkg = 0; iBkg < Nbkg; iBkg++ )
	{
	  TH1F* nomH = (TH1F*)dataFile.Get(bkgNames[iBkg]);
	  if( nomH == 0 )
	    {
	      std::cout << "cannot find histogram " << bkgNames[iBkg] << std::endl;
	      return;
	    }
	  bkgNominal.at(iBin).push_back(nomH->GetBinContent(iBin+1)+shift);
	  
	  if( iBkg == 0 ) {
	    scaleFactorsError.push_back( scaleFactors.at(iBin)* (SFErrHist->GetBinContent(iBin+1)/bkgNominal.at(iBin).at(0)-1) );
	    //std::cout << scaleFactors.at(iBin) << " +- " << scaleFactorsError.at(iBin) << std::endl;
	    float sfH = scaleFactors.at(iBin)/totalStatisticsHighLow.at(iBin).first*totalStatistics.at(iBin);
	    float sfL = scaleFactors.at(iBin)/totalStatisticsHighLow.at(iBin).second*totalStatistics.at(iBin);
	    
	    float sfHE = sfH/scaleFactors.at(iBin)*TMath::Sqrt( TMath::Power(scaleFactorsError.at(iBin),2)
								- 1/(totalStatistics.at(iBin))
								+ 1/totalStatisticsHighLow.at(iBin).first 
								);
	    float sfHE_v2 = sfH*TMath::Sqrt( TMath::Power(scaleFactorsError.at(iBin)/scaleFactors.at(iBin),2)
					     + 1/(totalStatistics.at(iBin))
					     + 1/totalStatisticsHighLow.at(iBin).first
					     );

	    float sfLE = sfL/scaleFactors.at(iBin)*TMath::Sqrt( TMath::Power(scaleFactorsError.at(iBin),2)
								- 1/(totalStatistics.at(iBin))
								+ 1/totalStatisticsHighLow.at(iBin).second);
	    std::cout << "iBin: " << iBin << std::endl;
	    std::cout << sfH << " +- " << sfHE << std::endl;
	    std::cout << sfH << " +- " << sfHE_v2 << std::endl;
	    std::cout << sfL << " +- " << sfLE << std::endl;
	    std::cout << std::endl;
	    
	    scaleFactorHigh.push_back(std::make_pair(sfH,sfHE));
	    scaleFactorLow.push_back(std::make_pair(sfL,sfLE));
	  }
	  //-----------------------------------------------------
	  //not clear what this does, besides add some systematic
	  //-----------------------------------------------------
	  //per systematic                             
	  for(int iSyst=0; iSyst<keys->GetEntries(); iSyst++) 
	    {
	      const char* name = keys->At(iSyst)->GetName();
	      //no match or its just the nominal
	      if( strncmp(bkgNames[iBkg].Data(), name,bkgNames[iBkg].Length()) != 0
		  || bkgNames[iBkg].Length() == strlen(name) ) 
		{
		  //delete name;                                
		  continue;
		}
	      TString tsName(name);
	      //delete name;
	      //only run for the Up systematic, we'll process both together to form the pair                                                          
	      if( !tsName.EndsWith("Up") ) 
		{
		  continue;
		}
	      
	      systMap systematics;
	      
	      TH1F* up   = (TH1F*)dataFile.Get(tsName);
	      TH1F* down = (TH1F*)dataFile.Get(tsName(0,tsName.Length()-2)+"Down");
	      if( up == 0 ) 
		{
		  std::cout << "could not find histogram named " << tsName << std::endl;
		  return;
		}
	      if( down == 0 )
		{
		  std::cout << "could not find histogram named " << tsName(0,tsName.Length()-2)+"Down" << std::endl;
		  return;
		}
	      
	      if( up->GetBinContent(iBin+1) == down->GetBinContent(iBin+1) ) continue;
	      //deal with the sideband statistcs systematic differently            
	      if( tsName.Contains("sidebandStatistics") || tsName.Contains("bkgShape") || tsName.Contains("bkg_fitUp") ) 
		{
		  continue;
		  float N;
		  //compute the statistics out of the card
		  if( up->GetBinContent(iBin+1)/bkgNominal.at(iBin).at(iBkg) == 1 ) N = 1;
		  else N = TMath::Power(1./(up->GetBinContent(iBin+1)/(bkgNominal.at(iBin).at(iBkg)-shift)-1),2);
		  float sf = (up->GetBinContent(iBin+1)+shift-bkgNominal.at(iBin).at(iBkg))/sqrt(N);
		  bkgStatistics.push_back( N );
		  //scaleFactors.push_back( sf);                        
		  std::cout << obsVec.at(iBin) << "  " << bkgNominal.at(iBin).at(iBkg) 
			    << "   " << up->GetBinContent(iBin+1) << "  " << N << "  " 
			    << sf << std::endl;
		  bkgSyst.at(iBin).at(iBkg).push_back( make_pair(bkgNominal.at(iBin).at(iBkg)+sqrt(N+0.76)*sf+shift,bkgNominal.at(iBin).at(iBkg)-sqrt(N	+0.76)*sf+shift) );
		  bkgSystNames.at(iBin).at(iBkg).push_back(tsName);
		}else
		{
		  bkgSyst.at(iBin).at(iBkg).push_back( make_pair(up->GetBinContent(iBin+1)+shift,down->GetBinContent(iBin+1)+shift));
		  bkgSystNames.at(iBin).at(iBkg).push_back(tsName);
		}
	    }
	}//end Bkg Loop
    }//end MR-Rsq bins loop

  TFile* _fout = new TFile("_profile_likelihood.root", "recreate");
  
  for( int iReg = 0; iReg < nBoxes; iReg++ )
    {
      vector<region> region = getSigRegions((BinningRegion)iReg);
      cout << endl <<endl << getRegionName( (BinningRegion)iReg ) << endl << endl;
      
      if(fullTex) 
	{
	  std::cout << "\\begin{table}\n\\begin{tabular}{|cc|c|c|cc|}" << std::endl
		    << "\\hline\n"
		    << "$M_R$ region & $R^2$ region & observed events & expected background & p-value & significance ($\\sigma$) \\\\\n" 
		    << "\\hline" << std::endl;
	}
      
    for( int iC = 0; iC < region.size(); iC++ ) 
      {
	int i = ConvertToCombineBin((BinningRegion)iReg,iC);
	float bkgTot = 0;
	for( int iBkg = 0; iBkg < Nbkg; iBkg++ ) bkgTot += bkgNominal.at(i).at(iBkg);
      
	float higgsTot = 0;
	pair<float,float> err = make_pair(0,0);
	for( int iBkg = 1; iBkg < Nbkg; iBkg++ ) //avoids non-resonant
	  { 
	    higgsTot += bkgNominal.at(i).at(iBkg);
	    pair<float,float> thisErr = addSystsQuad(bkgNominal.at(i).at(iBkg),bkgSyst.at(i).at(iBkg));
	    err.first  += thisErr.first;
	    err.second += thisErr.second;
	  }
	//----------------------------------------
	//total SM higgs bkg and its uncertainty
	//----------------------------------------
	std::pair<float,float> higgs = make_pair(higgsTot,sqrt(err.first));
	
	//----------------------------------------
	//nominal sideband non-resonant prediction
	//----------------------------------------
	double non_res_nominal_pred = scaleFactors.at(i)*bkgStatistics.at(i);
	//---------------------------------------------------------------------
	// non-resonant prediction statistical uncertainty = sf*sqrt(Nsideband)
	//---------------------------------------------------------------------
	double non_res_stat_err = scaleFactors.at(i)*TMath::Sqrt( bkgStatistics.at(i) );
	
	//----------------------------------------------------------
	//uppper and lower sideband non-resonant predictions
	//----------------------------------------------------------
	double non_res_upper_pred = scaleFactorHigh.at(i).first*bkgStatisticsHighLow.at(i).first;
	double non_res_lower_pred = scaleFactorLow.at(i).first*bkgStatisticsHighLow.at(i).second; 
	
	//---------------------------------------------------------------------------------------------
	//estimate the size of the deviation between uppper and lower sideband non-resonant predictions
	//---------------------------------------------------------------------------------------------
	double sideband_deviation_size = fabs( non_res_upper_pred - non_res_lower_pred );


	//-----------------------------------------------------------
	//compare size of the stat. uncertainty to sideband_deviation
	//extra_err is the size of an extra shape sytematic uncer.
	//-----------------------------------------------------------
	double extra_err = 0.0;
	if ( sideband_deviation_size > non_res_stat_err )
	  {
	    double plus_sigma  = fabs( max( non_res_upper_pred, non_res_lower_pred ) - non_res_nominal_pred );
	    double minus_sigma = fabs( non_res_nominal_pred - min( non_res_upper_pred, non_res_lower_pred ) );
	    extra_err = max( plus_sigma, minus_sigma );
	  }else
	  {
	    extra_err = non_res_stat_err;
	  }
	
	if( non_res_nominal_pred != 0 )
	  {
	    extra_err = 0.5*extra_err/non_res_nominal_pred;
	  }
	else
	  {
	    extra_err = 0.01;
	  }

	if ( extra_err > 0.5 ) extra_err = 0.5;//fixing maximum extra systematic error to 50%
	
	//-----------------------------------------
	//prepare input to minuit minimization
	//-----------------------------------------
	float obs = obsVec.at(i);
	std::pair<float,float> SF   =  make_pair(scaleFactors.at(i),scaleFactorsError.at(i));
	float n_sideband = bkgStatistics.at(i);
	extra_err = 0.5*extra_err/(float)n_sideband;

	//----------------------------------------
	//getting significance
	//----------------------------------------
	float sig = getSigMinuit(obs, n_sideband, SF, higgs, extra_err);
	if(isinf(sig)) return;
	float pv = RooStats::SignificanceToPValue(sig);

	//----------------------------------------                                            
	//getting delta log likelihood (for signal)                                                               
        //---------------------------------------- 
	TH1D* _dll_tmp = getDeltaLogLikelihood(obs, n_sideband, SF, higgs, extra_err, false);
	TString _h_name = Form("delta_log_likelihood_%d", i);
	_dll_tmp->Write( _h_name );
	
	
	//----------------------------------------
	//getting delta log likelihood (for Nobs)
	//----------------------------------------                                                                        
	TH1D* _dll_tmp_obs = new TH1D( *getDeltaLogLikelihood(obs, n_sideband, SF, higgs, extra_err, true) );
        _h_name = Form("delta_log_likelihood_obs_%d", i);
	std::pair<float, float> bkg_total_err = findOneSigma( _dll_tmp_obs );
	_dll_tmp_obs->Write( _h_name );
	
	//std::cout << "====> iBin: " << i << std::endl;
	//std::cout << "extra err = " << extra_err << std::endl;
	printf( "% 6.0f - % 6.0f & %0.2f - %0.2f & % 4.0f & $% 4.1f^{+%0.2f}_{-%0.2f}$ & %0.3f & %0.1f \\\\\n",
		region.at(iC).MR_min, region.at(iC).MR_max,
                region.at(iC).Rsq_min, region.at(iC).Rsq_max,
                obsVec.at(i), bkgTot, bkg_total_err.second, bkg_total_err.first,
		pv, sig
		);
	
	/*
	std::cout << "====> iBin: " << i << std::endl;
	std::cout << "Nobs: " << obsVec.at(i) << " Nexp: " << bkgTot << " (+" << bkg_total_err.second << ", -" << bkg_total_err.first 
		  << ")"
		  << " nH: " << higgs.first << " +/- " << higgs.second 
		  << " SF: " << scaleFactors.at(i) << " SFerr: " << scaleFactorsError.at(i)
		  << std::endl;
	std::cout << "nominal pred: " << non_res_nominal_pred
	<< ", upper-sideband pred: " << non_res_upper_pred
		  << ", lower-sideband pred: " << non_res_lower_pred
		  << std::endl;
	std::cout << "stat. uncertainty: " << non_res_stat_err
		  << ", |sidebandUpper-sidebandLower| =  " << sideband_deviation_size
		  << "; Extra uncertainty: " << extra_err
		  << std::endl;

	std::cout << "nsigmas: " << sig << ", p-val: " << pv 
	<< std::endl;
	*/
      }
    
    if( fullTex )
      {
	std::cout << "\\hline\n\\end{tabular}" << std::endl
		  << "\\caption{Number of Events observed in the signal region compared to expected"
		  << " background in the {\\bf " << getRegionName( (BinningRegion)iReg ) << "} box.}\n"
		  << "\\label{tab:ObsPV_" << getRegionName( (BinningRegion)iReg ) << "}\n"
		  << "\\end{table}" << std::endl;
      }
    
    }
  
  _fout->Close();
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

float getSigMinuit(float obs, float sideband, std::pair<float,float> SF, std::pair<float,float> higgsBkg,float bkgShapeErr)
{
  params.obs = obs;
  params.Nside = sideband;
  params.sf = SF.first;
  params.sfe = SF.second;
  params.Higgs = higgsBkg.first;
  params.HiggsErr = higgsBkg.second;
  params.combAddErr = bkgShapeErr;
  std::pair<float,float> prof = profileSimpleNoHistMinuit();
  
  return sqrt(prof.first);
};

TH1D* getDeltaLogLikelihood(float obs, float sideband, std::pair<float,float> SF,
			    std::pair<float,float> higgsBkg,float bkgShapeErr, bool _profileNobs )
{
  params.obs = obs;
  params.Nside = sideband;
  params.sf = SF.first;
  params.sfe = SF.second;
  params.Higgs = higgsBkg.first;
  params.HiggsErr = higgsBkg.second;
  params.combAddErr = bkgShapeErr;
  TH1D* _dll;
  if ( !_profileNobs )
    {
      _dll = getTwoLogLikelihood();
    }
  else
    {
      _dll = getTwoLogLikelihood( true );
    }

  return _dll;
};

pair<float, float> findOneSigma( TH1D* _nll )
{
  if ( _nll == NULL )return make_pair( 0, 0);
  float pSigma = 0.0, mSigma = 0.0;
  int nbins = _nll->GetNbinsX();
  float x_low = _nll->GetXaxis()->GetXmin();
  float step = ( _nll->GetXaxis()->GetXmax() - _nll->GetXaxis()->GetXmin() )/(float)nbins;
  int min_bin = _nll->GetMinimumBin();
  float x_min = x_low + step*(float)min_bin;
  float min_diff = 99999;
  for ( int i = 1; i <= min_bin; i++ )
    {
      if ( fabs( _nll->GetBinContent(i) - 1.0 ) < min_diff )
	{
	  mSigma = x_min - (x_low + step*(float)i);
	  min_diff = fabs( _nll->GetBinContent(i) - 1.0 );
	}
    }

  min_diff = 99999;
  for ( int i = min_bin; i <= nbins; i++ )
    {
      if ( fabs( _nll->GetBinContent(i) - 1.0 ) < min_diff )
        {
          pSigma = (x_low + step*(float)i) - x_min;
          min_diff = fabs( _nll->GetBinContent(i) - 1.0 );
        }
    }
  
  return make_pair( mSigma, pSigma);
};
