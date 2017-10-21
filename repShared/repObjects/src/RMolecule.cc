/// RMolecule.cc
/// ----------------------------------------------------------------------------
/// rep RMolecule c++ class 
/// Pierre GRANDEMANGE
/// 12/06/2017
/// ----------------------------------------------------------------------------



/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <iomanip>    /// fmtflags...
#include <cmath>
#include <stdlib.h>   /// getenv
#include <string>
#include <vector>
/// root classes
#include "TGraphErrors.h"
#include "TFile.h"
//#include "TCanvas.h"
#include "TPad.h"
/// rep namespaces
#include "RDump.h"
#include "RUnits.h"
/// rep classes
#include "RAtom.h"
#include "RFlyer.h"     /// needed for: ruth scat cross section methods
#include "RMachine.h"   /// needed for: ruth scat cross section methods
#include "RMolecule.h" 
#include "RParticle.h"
//#include "Machine.h"
using namespace std;


/// static attributes
/// **************************************************************************
/// fpbar ionisation cs file name
string RMolecule::rdataFolder( 
  (string)getenv("repData") + (string)"/crossSections/rootData/" );  
string RMolecule::rionCSfName( "tfPbarCrossSections.root" );
  


/// constructors, destructor, copy 
/// **************************************************************************
RMolecule::RMolecule(string const& name, string const& symbol)
  : RParticle(name, symbol), rcomposition(0), rsymbolIsotope("")
{}

RMolecule::RMolecule(string const& name, string const& symbol, 
                     string const& symbolWithIsotopes)
  : RParticle(name, symbol), rcomposition(0), rsymbolIsotope(symbolWithIsotopes)
{}

RMolecule::RMolecule(string const& name, string const& symbol, 
                     string const& symbolWithIsotopes, 
                     vector<RAtom*> const& composition)
  : RParticle(name, symbol), rcomposition(composition), 
    rsymbolIsotope(symbolWithIsotopes)
{}

RMolecule::~RMolecule()
{}



/// accessors
/// **************************************************************************
//vector<RAtom*> RMolecule::getComposition() const { return rcomposition; }



/// methods
/// ****************************************************************************
  
/// ion_pbarcs
/// ----------------------------------------------------------------------------
/// evaluate pbar cross section for given flyer(ekin)
/// using ion_pbarcs_data in this class
///   -> original graph units: 1.e-20m2 ; keV
///   -> input units: RUnits system (SI) -> default=Joule
///   -> outut units: RUnits system (SI) -> default=m2
///   -> interpolation is linear, ROOT TGraph "eval" method
double RMolecule::ion_pbarcs(RFlyer* const& flyer) const
{
  TGraphErrors* g=ion_pbarcs_data(); /// get ionsition cross section graph
  double eval = g->Eval(flyer->ekin()/RUnits::keV);
  return eval*1.e-20*RUnits::meter;
}
double RMolecule::ion_pbarcs(double ekin) const
{
  TGraphErrors* g=ion_pbarcs_data(); /// get ionsition cross section graph
  double eval = g->Eval(ekin/RUnits::keV);
  return eval*1.e-20*RUnits::meter;
}

/// ion_pbarcs_data
/// ----------------------------------------------------------------------------
/// retrieve TGraphErrors with pbar ionisation cross section(ekin) data
/// data from : "J. Phys. B: At. Mol. Opt. Phys. 44 (2011) 122001"
/// data stored and retrieved in a root file (also in .txt file)
/// units:  kinE [kEV]
///         cs   [1.e-20 m^2]
TGraphErrors* RMolecule::ion_pbarcs_data() const
{
  TGraphErrors* g(0);  // OUTPUT
    
  /// open TFile with cross-sections  
  const string tfPath = rdataFolder + rionCSfName;  /// static attibutes
  TFile tf(tfPath.c_str(), "READ");       /// open root file
  if( tf.IsZombie() ) { /// test if file is open
    cout << endl;
    cout << "------------------------------------- " << endl;
    cout << "rep warning: error opening the TFile: " << endl;
    cout << tfPath << endl;
    cout << "------------------------------------- " << endl;
    g = new TGraphErrors("emptyGraph");
  }
  else {
    //cout << "found the TFile" << endl; 
    string gName=rsymbol + "--All";
    g = (TGraphErrors*)tf.Get(gName.c_str());
  }
  tf.Close(); /// check if pointer graph pointer still valid when closing file
  return g;
}







/// rs_ebucs
/// ----------------------------------------------------------------------------
/// rutherford scattering total cross section approximation
/// IPAC2014, C.Carli, TUPRI028
/// sum of atomic cross sections
double RMolecule::rs_ebucs(RFlyer* const& flyer, RMachine* const& machine) const
{
  double dbl(0);
  for( int i=0; i<rcomposition.size(); ++i ) { 
    dbl += rcomposition[i]->rs_ebucs(flyer, machine);
  }
  return dbl;
}

/// rs_lalcs
/// ----------------------------------------------------------------------------
/// rutherford scattering large angle loss cross section
/// IPAC2014, C.Carli, TUPRI028
/// sum of atomic cross sections
double RMolecule::rs_lalcs(RFlyer* const& flyer, RMachine* const& machine) const
{
  double dbl(0);
  for( int i=0; i<rcomposition.size(); ++i ) { 
    dbl += rcomposition[i]->rs_lalcs(flyer, machine);
  }
  return dbl;
}

/// rs_totcs
/// ----------------------------------------------------------------------------
/// rutherford scattering total cross section approximation
/// IPAC2014, C.Carli, TUPRI028.
/// sum of atomic cross sections
double RMolecule::rs_totcs(string const& model) const
{
  double csMol(0);
  for( int i=0; i<rcomposition.size(); ++i ) { 
    csMol += rcomposition[i]->rs_totcs(model);
  }
  return csMol;
}





/// pKfactor
/// --------------------------------------------------------------------------
/// get pressure correction factor from pN2eq corresponding to gauge model 
/// data from pfeiffer pkr261 gauge manual
double RMolecule::pKfactor(double pN2eq, string const& gauge, 
													 bool verbose) const
{
  double pkf(1.);
  if( gauge=="pkr261" && pN2eq/RUnits::Pa < 1. ) {
    if      ( this->getSymbol()=="N2" )   pkf = 1.;
    else if ( this->getSymbol()=="O2" )   pkf = 1.;
    else if ( this->getSymbol()=="CO" )   pkf = 1.;
    else if ( this->getSymbol()=="Xe" )   pkf = 0.4;
    else if ( this->getSymbol()=="Kr" )   pkf = 0.5;
    else if ( this->getSymbol()=="Ar" )   pkf = 0.8;
    else if ( this->getSymbol()=="H2" )   pkf = 2.4;
    else if ( this->getSymbol()=="Ne" )   pkf = 4.1;
    else if ( this->getSymbol()=="He" )   pkf = 5.9;
    else {
			if( verbose==true ) {
				cout << "rep WARNING: no pK correction factor for molecule: ";
				cout << this->getSymbol() << endl;
			}	
    }
  } 
  else {
  	if( verbose==true ) {
  		cout << "rep WARNING: no pK correction factor for this p/gauge:";
  		cout << pN2eq << "/" << gauge << endl;
  	}
  }
  return pkf;
}



/// dump properties
/// **************************************************************************
/// dump
/// --------------------------------------------------------------------------
/// dump properties
void RMolecule::dump(ostream &flux) const
{
  flux << "RMolecule dump method" << endl; 
  RDump::line(21, "-", flux); 
  flux << rname << "  " << rsymbol << "   ";
  flux << rsymbolIsotope << "  " << rmass/RUnits::amu << "amu" << "  ";
  flux << rcharge/RUnits::e << "e" << "  ";
  flux << endl;
  RDump::line(46, "-", flux); 
  for( int i=0; i<rcomposition.size(); ++i) { rcomposition[i]->dumpLine(flux); }
  
}

/// dumpLine
/// --------------------------------------------------------------------------
  void RMolecule::dumpLine(ostream &flux) const
{
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(3); flux << scientific; flux << left;
  
  flux << rname             << "  ";
  flux << rsymbol           << "  ";
  flux << rsymbolIsotope    << "  ";
  flux << rmass/RUnits::amu << "amu"  << "  ";
  flux.precision(1); flux << fixed;
  flux << rcharge/RUnits::e << "e"    << "  ";
  flux << endl;
  RDump::line(51, "-", flux); 
  
  flux.flags(f);  /// restore flags
  
  /// dump atoms in molecules
  for( int i=0; i<rcomposition.size(); ++i) { 
    flux << "  -> "; rcomposition[i]->dumpLine(flux); 
  }
}

/// dumpcs
/// --------------------------------------------------------------------------
void RMolecule::dumpcs(RFlyer* fly, RMachine* mac, std::ostream &flux) const
{
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  int ll(80);
  
  flux << "RMolecule dumpcs method ";
  flux << "(cross sections: ionisation and rutherford scattering)" << endl; 
  RDump::line(ll, "-", flux); 
  flux << "mol: "; this->dumpLine();
  flux << "flyer: "; fly->dumpLine();
  flux << "       "; fly->dumpBeta();
  flux << "machine: "; mac->dumpLine();
  RDump::line(ll, "-", flux); 
  
  flux.precision(2); flux << scientific;
  flux << "ion_pbarcs = " << ion_pbarcs(fly)/RUnits::m2 << "m2" << " | ";
  flux << "ionisation";
  flux << endl;
  
  flux.precision(2); flux << scientific;
  flux << "rs_totcs =   " << rs_totcs("carli2014")/RUnits::m2 << "m2" << " | ";
  flux << "rutherford scattering - total";
  flux << endl;
  
  flux << "rs_lalcs =   " << rs_lalcs(fly,mac)/RUnits::m2 << "m2" << " | ";
  flux << "rutherford scattering - large angle loss";
  flux << endl;
  
  flux << "rs_ebucs =   " << rs_ebucs(fly,mac)/RUnits::m2 << "m2" << " | ";
  flux << "rutherford scattering - emittance blowup";
  flux << endl;
  
  RDump::line(ll, "-", flux); flux.flags(f);  /// restore flags
  
  /*
  flux << "molecule cross sections : " << this->getSymbol() << endl;
  for(int i=0; i<35; ++i) {flux << "-";} flux << endl; /// line
  flux << "total scat cs (carli2014): " << rs_totcs("carli2014") << endl;
  flux << "large angle scat loss    : " << rs_lalcs(fly,mac) << endl;
  flux << "emittance blow-up cs     : " << rs_ebucs(fly,mac) << endl;
  flux << "pbar ionisation cs       : " << ion_pbarcs(fly) << endl;
  for(int i=0; i<35; ++i) {flux << "-";} flux << endl; /// line
  */
}


