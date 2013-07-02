#include "include/ElectronChecks.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "CommonTools/include/LeptonIdBits.h"

ElectronChecks::ElectronChecks(TTree *tree) : Vecbos(tree) 
{
  // common kinematic selections
  std::string theConfigDir       = "config/vecbos/";
  std::string fileCutsCommon     = theConfigDir + "ElectronChecksCuts.txt";
  std::string fileSwitchesCommon = theConfigDir + "ElectronChecksSwitches.txt";

  _commonSel = new Selection(fileCutsCommon,fileSwitchesCommon);
  ConfigCommonSelections(_commonSel);

  std::cout << "[GoodRunLS]::goodRunLS is " << _commonSel->getSwitch("goodRunLS") << " isData is " <<  _commonSel->getSwitch("isData") << std::endl;

  //To read good run list!
  if (_commonSel->getSwitch("goodRunLS") && _commonSel->getSwitch("isData"))
    {
      std::string goodRunGiasoneFile       = "config/vecbos/json/goodRunLS.json";
      setJsonGoodRunList(goodRunGiasoneFile);
      fillRunLSMap();
    }

  cutLevel[0]="raw";
  cutLevel[1]="rawSpikeId";
  cutLevel[2]="etGt5SpikeId";
  cutLevel[3]="etGt5LooseId";

  category[0]="all";
  category[1]="all_eb";
  category[2]="all_ee";
  category[3]="ed";
  category[4]="ed_eb";
  category[5]="ed_ee";
  category[6]="td";
  category[7]="td_eb";
  category[8]="td_ee";

  bookHistos();
}

ElectronChecks::~ElectronChecks(){}

void ElectronChecks::fillElectronPlots(unsigned int iele, unsigned int icutLevel, unsigned int icat)
{
  h_ele_pt[icutLevel][icat]->Fill(etEle[iele]);
  h_ele_eta[icutLevel][icat]->Fill(etaEle[iele]);
  h_ele_phi[icutLevel][icat]->Fill(phiEle[iele]);

  //Filling eleId variables
  float HoE, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[iele], bits::isEcalDriven);
  HoE = hOverEEle[iele];
  deta = deltaEtaAtVtxEle[iele];
  dphiin = deltaPhiAtVtxEle[iele];
  dphiout = deltaPhiAtCaloEle[iele];
  fbrem = FBrem(iele);
  eopout = eSeedOverPoutEle[iele];
  eop = eSuperClusterOverPEle[iele];
  if(ecaldriven) 
    {
      int sc = superClusterIndexEle[iele];
      see = sqrt(covIEtaIEtaSC[sc]);
      spp = sqrt(covIPhiIPhiSC[sc]);
    } 
  else 
    {
      int sc = PFsuperClusterIndexEle[iele];
      if(sc>-1) 
	{
	  see = sqrt(covIEtaIEtaPFSC[sc]);
	  spp = sqrt(covIPhiIPhiPFSC[sc]);
	} 
      else 
	{
	  see = 999.;
	  spp = 999.;
	}
    }  

  h_ele_hOverE[icutLevel][icat]->Fill(HoE);
  h_ele_sigmaIetaIeta[icutLevel][icat]->Fill(see);
  h_ele_sigmaIphiIphi[icutLevel][icat]->Fill(spp);
  h_ele_dEtaIn[icutLevel][icat]->Fill(deta);
  h_ele_dPhiIn[icutLevel][icat]->Fill(dphiin);
  h_ele_dPhiOut[icutLevel][icat]->Fill(dphiout);
  h_ele_fBrem[icutLevel][icat]->Fill(fbrem);
  h_ele_eOverP[icutLevel][icat]->Fill(eop);
  h_ele_eSeedOverPOut[icutLevel][icat]->Fill(eopout);
  h_ele_class[icutLevel][icat]->Fill(classificationEle[iele]+1.);

  float EWKz = trackVzGsfTrack[iele];
  
  // search for PV
  int closestPV = -1;   
  float z0 = 0.0;
  if(nPV>0) 
    {         // a primary vertex was found
      float minDzPV=999.;
      for(int v=0; v<nPV; v++) 
	{
	  if(fabs(PVzPV[v]-EWKz)<minDzPV) 
	    {
	      minDzPV=fabs(PVzPV[v]-EWKz);
	      closestPV = v;
	    }
	}
      z0 = PVzPV[closestPV];    // if nPV=0 --> z0 = 0
    }
  
  int gsfTrack = gsfTrackIndexEle[iele];
  // zele - zPV
  float dxy,dz;

  if (closestPV>-1)
    {
      dxy = eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], trackVxGsfTrack[gsfTrack], trackVyGsfTrack[gsfTrack], trackVzGsfTrack[gsfTrack], pxEle[iele], pyEle[iele], pzEle[iele]);
      dz = trackVzGsfTrack[gsfTrack]-PVzPV[closestPV];
    }
  else
    {
      dxy = eleDxyPV( 0., 0. , 0. , trackVxGsfTrack[gsfTrack], trackVyGsfTrack[gsfTrack], trackVzGsfTrack[gsfTrack], pxEle[iele], pyEle[iele], pzEle[iele]);
      dz = trackVzGsfTrack[gsfTrack]-0.;
    }
  
  //Tip Lip
  h_ele_tip[icutLevel][icat]->Fill(dxy);
  h_ele_lip[icutLevel][icat]->Fill(dz);

  //Isolation
  h_ele_dr04EJIso[icutLevel][icat]->Fill( dr04EcalRecHitSumEtEle[iele] );
  h_ele_dr04HIso[icutLevel][icat]->Fill( dr04HcalTowerSumEtEle[iele] );
  h_ele_dr04TrIso[icutLevel][icat]->Fill( dr04TkSumPtEle[iele] ) ;
}

bool ElectronChecks::isEleSelected(unsigned int iele,const TString& cutLevel)
{
  if (cutLevel == "raw" )
    {
      return true;
    }
  else if (cutLevel == "rawSpikeId" )
    {
      float rawEtSC=rawEnergySC[superClusterIndexEle[iele]]/TMath::CosH(etaSC[superClusterIndexEle[iele]]);
      if ( rawEtSC > 4. && 
	   eMaxSC[superClusterIndexEle[iele]]/e3x3SC[superClusterIndexEle[iele]]> 0.95 
	   )
	return false;

      return true;
    }
  else if (cutLevel == "etGt10SpikeId" )
    {
      if ( etEle[iele] < 10. ) 
	return false;
      if (  eMaxSC[superClusterIndexEle[iele]]/e3x3SC[superClusterIndexEle[iele]]> 0.95 )
	return false;
      
      return true;
    }
  else if (cutLevel == "etGt10LooseId" )
      {
	if ( etEle[iele] < 10. ) 
	  return false;
	if (  eMaxSC[superClusterIndexEle[iele]]/e3x3SC[superClusterIndexEle[iele]]> 0.95 )
	  return false;
	if ( ! (eleIdCutsEle[iele]& 0x4) )
	  return false;
	
	return true;
      }
  
  
};

bool ElectronChecks::isEleCategory(unsigned int iele,unsigned int icat)
{

  bool isEB = fiducialFlagsEle[iele] & 0x200;
  bool isEE = fiducialFlagsEle[iele] & 0x100;
  bool isED = recoFlagsEle[iele] & 0x2;

  //All
  if (icat == 0)
    {
      return true;
    }
  //EB
  else if (icat == 1)
    {
      if (isEB)
	return true;
    }
  //EE
  else if (icat == 2)
    {
      if (isEE)
	return true;
    }
  //Ecal Driven
  else if (icat == 3)
    {
      if (isED)
	return true;
    }
  //Ecal Driven EB
  else if (icat == 4)
    {
      if ( isED && isEB)
	return true;
    }
  //Ecal Driven EE
  else if (icat == 5)
    {
      if ( isED && isEE)
	return true;
    }
  //Tracker Driven Only
  else if (icat == 6)
    {
      if ( !isED )
	return true;
    }
  //Tracker Driven Only EB
  else if (icat == 7)
    {
      if ( !isED && isEB )
	return true;
    }
  //Tracker Driven Only EE
  else if (icat == 8)
    {
      if ( !isED && isEE )
	return true;
    }

  return false;
}

bool ElectronChecks::isEventSelected()
{
  return true;
}

void ElectronChecks::Loop()
{

   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "Total number of entries in the chain = " << nentries << std::endl;   
   Long64_t nbytes = 0, nb = 0;
   
   unsigned int lastLumi=0;
   unsigned int lastRun=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);

     if (ientry < 0) break;
     if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

     //Good runLS
     if ( _commonSel->getSwitch("goodRunLS") && _commonSel->getSwitch("isData") && !isGoodRunLS())
       {
	 if ( lastRun!= runNumber || lastLumi != lumiBlock)
	   {
	     lastRun = runNumber;
	     lastLumi = lumiBlock;
	     std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
	   }
	 continue;
       }
     
     if ( _commonSel->getSwitch("goodRunLS") && _commonSel->getSwitch("isData") && ( lastRun!= runNumber || lastLumi != lumiBlock) )
       {
	 lastRun = runNumber;
	 lastLumi = lumiBlock;
	 std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
       }
     
     if (!isEventSelected())
       continue;
     
     //Analysis
     for (unsigned int iele=0;iele<nEle;++iele)
       {
	 for (unsigned int icut=0;icut<NCUTLEVELS;++icut)
	   {
	     for (unsigned int icat=0;icat<NCATEGORIES;++icat)
	       {
		 if (!isEleSelected(iele,cutLevel[icut]))
		   continue;
		 if (!isEleCategory(iele,icat))
		   continue;
		 fillElectronPlots(iele,icut,icat);
	       }
	   }
       }
   }
   
}

void ElectronChecks::ConfigCommonSelections(Selection* _selection) 
{
  _selection->addSwitch("isData");
  _selection->addSwitch("goodRunLS");
  _selection->summary();
}

void ElectronChecks::bookHistos()
{
   for (unsigned int icut=0;icut<NCUTLEVELS;++icut)
     {
       for (unsigned int icat=0;icat<NCATEGORIES;++icat)
	 {
	   //Kinematics
	   h_ele_pt[icut][icat]=book1D("h_ele_pt_"+cutLevel[icut]+"_"+category[icat],"h_ele_pt_"+cutLevel[icut]+"_"+category[icat],200,0.,200.);
	   h_ele_eta[icut][icat]=book1D("h_ele_eta_"+cutLevel[icut]+"_"+category[icat],"h_ele_eta_"+cutLevel[icut]+"_"+category[icat],60,-3.,3.);
	   h_ele_phi[icut][icat]=book1D("h_ele_phi_"+cutLevel[icut]+"_"+category[icat],"h_ele_phi_"+cutLevel[icut]+"_"+category[icat],60,-TMath::Pi(),TMath::Pi());

	   //EleId Variables
	   h_ele_hOverE[icut][icat]=book1D("h_ele_hOverE_"+cutLevel[icut]+"_"+category[icat],"h_ele_hOverE_"+cutLevel[icut]+"_"+category[icat],50,0.,0.3);
	   h_ele_sigmaIetaIeta[icut][icat]=book1D("h_ele_sigmaIetaIeta_"+cutLevel[icut]+"_"+category[icat],"h_ele_sigmaIetaIeta_"+cutLevel[icut]+"_"+category[icat],200,0.,0.06);
	   h_ele_sigmaIphiIphi[icut][icat]=book1D("h_ele_sigmaIphiIphi_"+cutLevel[icut]+"_"+category[icat],"h_ele_sigmaIphiIphi_"+cutLevel[icut]+"_"+category[icat],200,0.,0.2);
	   h_ele_dEtaIn[icut][icat]=book1D("h_ele_dEtaIn_"+cutLevel[icut]+"_"+category[icat],"h_ele_dEtaIn_"+cutLevel[icut]+"_"+category[icat],100,-0.05,0.05);
	   h_ele_dPhiIn[icut][icat]=book1D("h_ele_dPhiIn_"+cutLevel[icut]+"_"+category[icat],"h_ele_dPhiIn_"+cutLevel[icut]+"_"+category[icat],100,-0.2,0.2);
	   h_ele_dPhiOut[icut][icat]=book1D("h_ele_dPhiOut_"+cutLevel[icut]+"_"+category[icat],"h_ele_dPhiOut_"+cutLevel[icut]+"_"+category[icat],100,-0.2,0.2);
	   h_ele_fBrem[icut][icat]=book1D("h_ele_fBrem_"+cutLevel[icut]+"_"+category[icat],"h_ele_fBrem_"+cutLevel[icut]+"_"+category[icat],110,-0.1,1.0);
	   h_ele_eOverP[icut][icat]=book1D("h_ele_eOverP_"+cutLevel[icut]+"_"+category[icat],"h_ele_eOverP_"+cutLevel[icut]+"_"+category[icat],250,0.,5.);
	   h_ele_eSeedOverPOut[icut][icat]=book1D("h_ele_eSeedOverPOut_"+cutLevel[icut]+"_"+category[icat],"h_ele_eSeedOverPOut_"+cutLevel[icut]+"_"+category[icat],250,0.,5.);
	   h_ele_class[icut][icat]=book1D("h_ele_class_"+cutLevel[icut]+"_"+category[icat],"h_ele_class_"+cutLevel[icut]+"_"+category[icat],6,-0.5,5.5);
	   //Tip Lip
	   h_ele_tip[icut][icat]=book1D("h_ele_tip_"+cutLevel[icut]+"_"+category[icat],"h_ele_tip_"+cutLevel[icut]+"_"+category[icat],200,-0.5,0.5);
	   h_ele_lip[icut][icat]=book1D("h_ele_lip_"+cutLevel[icut]+"_"+category[icat],"h_ele_lip_"+cutLevel[icut]+"_"+category[icat],200,-0.5,0.5);

	   //Isolation
	   h_ele_dr04EJIso[icut][icat]=book1D("h_ele_dr04EJIso_"+cutLevel[icut]+"_"+category[icat],"h_ele_dr04EJIso_"+cutLevel[icut]+"_"+category[icat],110,-1.,10.);
	   h_ele_dr04HIso[icut][icat]=book1D("h_ele_dr04HIso_"+cutLevel[icut]+"_"+category[icat],"h_ele_dr04HIso_"+cutLevel[icut]+"_"+category[icat],110,-1.,10.);
	   h_ele_dr04TrIso[icut][icat]=book1D("h_ele_dr04TrIso_"+cutLevel[icut]+"_"+category[icat],"h_ele_dr04TrIso_"+cutLevel[icut]+"_"+category[icat],110,-1.,10.);

	 }
     }
}
