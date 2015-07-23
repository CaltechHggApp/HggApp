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
  std::cout.precision(3);
    
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
      std::cout << "\\begin{table}\n\\begin{tabular}{|c|c|c|c|c|}\n\\hline\n";
      std::cout << "Box & $N_{1}$ & $N_{2}$ & $\\alpha_{1}$ & $\\alpha_{2}$\\\\\n\\hline";
    }
  for ( int i = 0; i < nCat; i++ )
    {
      RooWorkspace* w = (RooWorkspace*)dataFile->Get( cat[i] + "_mgg_workspace" );
      if( _fullTable )
	{
	  std::cout << "\n" << cat[i] << " & " 
		    << w->var("Nbkg1")->getVal()*w->var("Nbkg1")->getVal() << " $\\pm$ " 
		    << 2.0*w->var("Nbkg1")->getError() << " & "
		    << w->var("Nbkg2")->getVal()*w->var("Nbkg2")->getVal() << " $\\pm$ "
                    << 2.0*w->var("Nbkg2")->getError() << " & "
		    << w->var("a1")->getVal()*w->var("a1")->getVal() << " $\\pm$ "
                    << 2.0*w->var("a1")->getError() << " & "
		    << w->var("a2")->getVal()*w->var("a2")->getVal() << " $\\pm$ "
                    << 2.0*w->var("a2")->getError() << "\\\\";	  
	}
      else
	{
	  std::cout << cat[i] << "\nN1:\t " << w->var("Nbkg1")->getVal()*w->var("Nbkg1")->getVal() << " +/- "
                    << 2.0*w->var("Nbkg1")->getError()
                    << "\nN2:\t " << w->var("Nbkg2")->getVal()*w->var("Nbkg2")->getVal() << " +/- "
                    << 2.0*w->var("Nbkg2")->getError()
                    << "\nalpha1:\t " <<w->var("a1")->getVal()*w->var("a1")->getVal() << " +/- "
                    << 2.0*w->var("a1")->getError() 
                    << "\nalpha2:\t "<< w->var("a2")->getVal()*w->var("a2")->getVal() << " +/- "
                    << 2.0*w->var("a2")->getError() << "\n";
	}
    }

  //--------------------------------
  // f u l l T a b l e   e n d i n g
  //--------------------------------
  if( _fullTable )
    {
      std::cout << "\n\\hline\n\\end{tabular}\n\\\caption{";
      std::cout << "\nParameters of the fits to the $m_{\\gamma\\gamma}$ upper and lower sideband. The quoted uncertai\
nties account for the statistical uncertainty in the fit}";
      std::cout << "\n\\label{tab:fit_parameters}\n\\end{table}\n";
    }
  
  return;
}
