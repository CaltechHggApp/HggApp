// This is used to make the background fit to data
// run as a root macro right now, takes a list to the output from HggSelectorApp

#include "RooRealVar.h"         
#include "RooArgSet.h"          
#include "RooArgList.h"         
#include "RooDataSet.h"         
#include "RooExponential.h"     
#include "RooLandau.h"          
#include "RooPlot.h"            
#include "RooFit.h"             
#include "RooAddPdf.h"          
#include "RooGaussian.h"        
#include "RooCBShape.h"         
#include "RooFFTConvPdf.h"      
#include "RooDataHist.h"        
#include "RooHistPdf.h"         
#include "RooHistFunc.h"        
#include "RooMoment.h"          
#include "RooFitResult.h"       
#include "RooExtendPdf.h"       
#include "RooGenericPdf.h"      
#include "RooBreitWigner.h"     
#include "RooBifurGauss.h"      
#include "RooProdPdf.h"         
#include "RooCategory.h"        
#include "RooSimultaneous.h"    
#include "RooWorkspace.h"       
#include "RooConstVar.h"        
#include "RooChebychev.h"       
#include "RooBernstein.h"

using namespace RooFit; 

struct TreeData{
   Int_t           trigger;
   Float_t         mPair;
   Float_t         mPairRes;
   Float_t         mPairResWrongVtx;
   Float_t         diPhotonMVA;
   Float_t         Photon1MVA;
   Float_t         Photon2MVA;
   Int_t           diPhotonVtx;
   Float_t         PFmPair;
   Float_t         PFmPairRes;
   Float_t         PFmPairResWrongVtx;
   Float_t         PFdiPhotonMVA;
   Float_t         PFPhoton1MVA;
   Float_t         PFPhoton2MVA;
   Int_t           PFdiPhotonVtx;
   Float_t         energyPho1;
   Float_t         etaPho1;
   Float_t         pxPho1;
   Float_t         pyPho1;
   Float_t         pzPho1;
   Float_t         indexPho1;
   Float_t         energyPho2;
   Float_t         etaPho2;
   Float_t         pxPho2;
   Float_t         pyPho2;
   Float_t         pzPho2;
   Float_t         indexPho2;
   Float_t         energyPFPho1;
   Float_t         etaPFPho1;
   Float_t         pxPFPho1;
   Float_t         pyPFPho1;
   Float_t         pzPFPho1;
   Float_t         indexPFPho1;
   Float_t         energyPFPho2;
   Float_t         etaPFPho2;
   Float_t         pxPFPho2;
   Float_t         pyPFPho2;
   Float_t         pzPFPho2;
   Float_t         indexPFPho2;
};

void MakeBackgroundFit(){
  

}
