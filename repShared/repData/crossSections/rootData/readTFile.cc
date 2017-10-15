// Read the TFile with cross sections generated by csToRoot
// Pierre Grandemange
// 12/02/2015

// standard library
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
// root classes
#include "TFile.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TMath.h"

using namespace std;

// Prototypes ******************************************************************
void readTFile();
double selectCS(string const& gasName, double const& energyKin);
double computeCrossSection_H2(double const energyKin);
//void drawCSused(TGraphErros csGraph);

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// Main program interpreted or compiled on the fly by CINT (ROOT interpreter)
void mainRoot()
{
  cout << endl;
  cout << "******************" << endl;
  cout << "Main program start" << endl;
  cout << "******************" << endl;
  cout << endl;

  readTFile();
  cout << endl;
  double cs=selectCS("CO2", 2000);
  cout << endl;
  cout << cs << endl;
  
  
  cout << endl;
  cout << "****************" << endl;
  cout << "Main program end" << endl;
  cout << "****************" << endl;
  cout << endl;
}

// *****************************************************************************
// main program compiled by g++ compiler
#ifndef __CINT__
int main()
{
  mainRoot();
  return 0;
}
#endif // __CINT__
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// Functions *******************************************************************
void readTFile()
{
  string const tfName="tfPbarCrossSections.root";   // rootfile name
  
  TFile *tf = new TFile(tfName.c_str(), "READ");
  if (tf->IsZombie()) { cout << "Error opening : " << tfName << endl; }
  else 
  { 
    tf->Print();
    tf->ls();
  }
}

// *****************************************************************************
// Select cross section from gas name and particle energy
// Interpolate cross section
double selectCS(string const& gasName, double const& energyKin)
{
// output
  double CSselected(-11.0);   // OUTPUT, selected cross section

// Open TFile with all cross sections  
  string const tfName="tfPbarCrossSections.root";   // rootfile name
  TFile *tf = new TFile(tfName.c_str(), "READ");    // open root file
  
  if (tf->IsZombie())                               // test if file is open
  { cout << "Error opening : " << tfName << endl; }
  else 
  {
// Select cross section by gas name    
    if (gasName=="CO2")
    {
      // open TGraph with cross sections associted to gas name
      TGraphErrors *tge=(TGraphErrors*)tf->Get("CO2--CO2plus_Knudsen");
      int const nPoints=tge->GetN();    // nb points in graph
      double *xG=tge->GetX();           // x array of graph
      double const xMax=TMath::MaxElement(nPoints, xG); // x array max element
      double const xMin=TMath::MinElement(nPoints, xG); // x array min element
   
      cout << nPoints << endl;
      cout << xMax << endl;
      cout << xMin << endl;
         
// test if particle energy is in the energy range of cross section data
      if ((energyKin>xMin) && (energyKin<xMax))
      {
        cout << "inside" << endl;
        CSselected = tge->Eval(energyKin);
      }
      else
      {
        cout << "!!!!!!!!!!!!!!!" << endl;
        cout << "!!! WARNING !!!" << endl;
        cout << "!!!!!!!!!!!!!!!" << endl;
        cout << "particle energy is out of range of data from : " << endl;
        cout << tge->GetTitle() << endl;
        cout << "The cross section is computed for proton and dihydrogen"
        << "by aproximation formulae from : " << endl;
        cout << "B.Hochadel et al. / Nucl. Instr. and Meth. in Phys. Res. A "
        << "343 (1994) 401-414" << endl;
        
        CSselected=computeCrossSection_H2(energyKin);
      }
    }
    else
    {
      cout << "gas unknown" << endl;
    }
  }
  
  return CSselected;
}


// *****************************************************************************
double computeCrossSection_H2(double const energyKin)
// Cross section calculation for a given energy - ion/dihydrogen
// B.Hochadel et al. / Nucl. Instr. and Meth. in Phys. Res. A 343 (1994) 401-414
{
  double crossSection(-11.);        // [10^-20 m^2] OUTPUT cross section
  double charge(1.);                // [e] particle charge
  double mass(1.0072765);           // pbar mass in atomic mass unit
  double energyKinU=energyKin/mass; // energy per atomic mass unit
  
  crossSection = 2*pow(charge,2)*(2*1.e-1)/energyKinU;
  
  return crossSection;
}

/*
// *****************************************************************************
void drawCSused(TGraphErros csGraph)
// Draw  multigraph of cross sections dataset/formulae used to select 
// cross section 
{
  TMultiGraph *tmg=new TMultiGraph();
}
*/
