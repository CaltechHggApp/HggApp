//C++ INCLUDES
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TString.h>


void makeScaleFactorTable( TString dir = "./", bool _fullTable = false )
{
  TFile* dataFile = new TFile( dir + "/data.root" );
  if ( !dataFile->IsOpen() )
    {
      std::cerr << "[ERROR]: CANNOT OPEN: " << dir + "/data.root" << std::endl;
      return;
    }
  
  //--------------------------------------------
  // s e t   i o s t r e a m   p r e c i s i o n
  //--------------------------------------------
  std::cout << std::fixed;
  std::cout.precision(4);
    
  //--------------------------------------
  // D e f i n i n g   c a t e g o r i e s
  //--------------------------------------
  int nCat = 5;
  TString cat[] = {"HighPt","Hbb", "Zbb", "HighRes", "LowRes"}; 
  
  //--------------------------------
  // f u l l T a b l e   H e a d e r
  //--------------------------------
  if( _fullTable )
    {
      std::cout << "\\begin{table}\n\\begin{tabular}{|c|c|}\n\\hline\n";
      std::cout << "Box & Background Prediction Scale Factor\\\\\n\\hline";
    }
  for ( int i = 0; i < nCat; i++ )
    {
      RooWorkspace* w = (RooWorkspace*)dataFile->Get( cat[i] + "_mgg_workspace" );
      if( _fullTable )
	{
	  std::cout << "\n" << cat[i] << " & " 
		    << w->var("scaleFactor")->getVal() << " $\\pm$ " 
		    << w->var("scaleFactor")->getError() << "\\\\";
	}
      else
	{
	  std::cout << cat[i] + ":\t " << w->var("scaleFactor")->getVal() << " +/- "  
		    << w->var("scaleFactor")->getError() << std::endl;
	}
    }

  //--------------------------------
  // f u l l T a b l e   e n d i n g
  //--------------------------------
  if( _fullTable )
    {
      std::cout << "\n\\hline\n\\end{tabular}\n\\\caption{";
      std::cout << "\nScale factors derived from the fits to the $m_{\\gamma\\gamma}$ upper and lower sideband. The quoted uncertai\
nties account for the statistical uncertainty in the fit and the choice of the functional form.}";
      std::cout << "\n\\label{tab:scale_factors}\n\\end{table}\n";
    }
  
  return;
}
