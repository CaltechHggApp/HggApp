//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef ElectronChecks_h
#define ElectronChecks_h

#include "VecbosBase.hh"
#include "Vecbos.hh"

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Utils.hh"

#define NCUTLEVELS 4
#define NCATEGORIES 9

using namespace std;

class ElectronChecks : public Vecbos{
public:

  ElectronChecks(TTree *tree=0); /// Class Constructor
  virtual ~ElectronChecks();     /// Class Destructor
  void Loop();

  TH1F* book1D( const char* name, const char* title, int nbins, float xmin, float xmax)
  {
    TH1F* temp=new TH1F(name,title,nbins,xmin,xmax);
    histoList.push_back((TObject*)temp);
    return temp;
  }
  
  TH2F* book2D( const char* name, const char* title, int nbinsX, float xmin, float xmax, int nbinsY, float ymin, float ymax)
  {
    TH2F* temp=new TH2F(name,title,nbinsX,xmin,xmax,nbinsY,ymin,ymax);
    histoList.push_back((TObject*)temp);
    return temp;
  }

  void writeAllHistos(const TString& outFile)
  {
    std::cout << "Saving histograms " << std::endl;
    TFile* eleHisto=new TFile(outFile,"RECREATE");
    eleHisto->cd();
    for (std::vector<TObject*>::const_iterator histo=histoList.begin(); histo!=histoList.end(); ++histo)
	(*histo)->Write();

    eleHisto->Write();
    eleHisto->Close();
  }

private:

  void ConfigCommonSelections(Selection* _selection);
  void fillElectronPlots(unsigned int iele, unsigned int icutLevel, unsigned int icat);
  bool isEleSelected(unsigned int iele,const TString& cutLevel);
  bool isEleCategory(unsigned int iele,unsigned int icat);
  bool isEventSelected();
  void bookHistos();

  // [4] selection levels 
  // [9] categories (all,allEB, allEE, EcalDriven, EcalDrivenEB, EcalDrivenEE, !EcalDriven, !EcalDrivenEB, !EcalDrivenEE)

  TString cutLevel[NCUTLEVELS];
  TString category[NCATEGORIES];

  //Kinematics
  TH1F* h_ele_pt[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_eta[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_phi[NCUTLEVELS][NCATEGORIES];

  //EleId
  TH1F* h_ele_hOverE[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_sigmaIetaIeta[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_sigmaIphiIphi[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_dEtaIn[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_dPhiIn[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_dPhiOut[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_fBrem[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_eOverP[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_eSeedOverPOut[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_class[NCUTLEVELS][NCATEGORIES];

  //Tip Lip
  TH1F* h_ele_tip[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_lip[NCUTLEVELS][NCATEGORIES];
  
  //Isolation
  TH1F* h_ele_dr04EJIso[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_dr04HIso[NCUTLEVELS][NCATEGORIES];
  TH1F* h_ele_dr04TrIso[NCUTLEVELS][NCATEGORIES];

  Utils anaUtils;
  
  std::vector< TObject* > histoList;
  
  Selection *_commonSel;
};
#endif
