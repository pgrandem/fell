/// RAtom.cc
/// ----------------------------------------------------------------------------
/// rep RAtom c++ class 
/// Pierre GRANDEMANGE
/// 12/06/2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
/// root classes
#include "TMath.h"
/// rep namespaces
#include "RDump.h"
#include "RUnits.h"
/// rep classes
#include "RAtom.h"
#include "RFlyer.h"
#include "RMachine.h"
#include "RParticle.h"
using namespace std;


/// constructor, destructor, copy 
/// ****************************************************************************
RAtom::RAtom() : RParticle(), rA(0), rZ(0)
{}

RAtom::RAtom(string name, string symbol) : RParticle(name, symbol), rA(0), rZ(0)
{}

RAtom::RAtom(string name, string symbol, int Z, int A) 
  : RParticle(name, symbol, 0., A*RUnits::amu), rA(A), rZ(Z)
{}

RAtom::RAtom(string name, string symbol, int Z, int A, double mass) 
  : RParticle(name, symbol, 0., mass), rA(A), rZ(Z)
{}

RAtom::RAtom(string name, string symbol, int Z, int A, 
             double mass, double charge) 
  : RParticle(name, symbol, charge, mass), rA(A), rZ(Z)
{}

RAtom::~RAtom() 
{}



/// methods
/// ****************************************************************************
  
/// ____________________________________________________________________________
/// rs : rutherford scattering methods.
/// -> IPAC2014, TUPRI028, Carli modified by me
  
/// rs_ebucl()
/// ----------------------------------------------------------------------------
/// rutherford scattering emittance blow up coulomb logarithm
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rs_ebucl(RFlyer* const& flyer, RMachine* const& machine) const
{ return log(rs_totra("carli2014")/rs_lalip(flyer,machine)); }

/// rs_ebucs()
/// ----------------------------------------------------------------------------
/// rutherford scattering large angle loss cross section
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rs_ebucs(RFlyer* const& flyer, RMachine* const& machine) const
{ return pow(rs_ebuip(flyer, machine),2)*M_PI; }

/// rs_ebuip()
/// ----------------------------------------------------------------------------
/// rutherford scattering large angle loss cross section
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rs_ebuip(RFlyer* const& flyer, RMachine* const& machine) const
{ return rs_lalip(flyer,machine) * sqrt( rs_ebucl(flyer,machine)/2. ); }




/// rs_lalcs()
/// ----------------------------------------------------------------------------
/// rutherford scattering large angle loss cross section
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rs_lalcs(RFlyer* const& flyer, RMachine* const& machine) const
{ return pow(rs_lalip(flyer, machine),2)*M_PI; }

/// rs_lalip()
/// ----------------------------------------------------------------------------
/// rutherford scattering large angle loss impact parameter
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rs_lalip(RFlyer* const& flyer, RMachine* const& machine) const
{
  double dbl;
  double Z(rZ);                             /// gas atom Z
  double z(flyer->getCharge()/RUnits::e);   /// flyer charge state, +/-1
  double m(flyer->getMass());               /// flyer mass
  double v(flyer->velocity());              /// flyer velocity
  double phi(machine->acceptanceAngle());   /// machine acceptance angle
  //cout << Z << "   " << z << "   " << m << "   " << m/RUnits::amu << "   ";
  //cout << v << "   " << phi << endl;
  dbl=rsip(Z, z, m, v, phi);  /// rutherford scattering formula
  return dbl;
}



/// rs_totcs()
/// ----------------------------------------------------------------------------
/// rutherford scattering total cross section approximation
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rs_totcs(string const& model) const
{
  double cs(0.);          /// [ ] cross section, OUTPUT
  if( model=="carli2014" ) {
    /// Procedings of IPAC2014, C.Carli, CORRECTED (e-21 in paper -> e-20 )
    cs = 0.79*1.e-20*pow(rZ,2./3.)*RUnits::m2;
  }
  /// Sternglass Phys.Rev.108-1
  //cs = 1.6 * pow(m_charge, 1./3.)*1.e-20*unitF;
  
  return cs;
}

/// rs_totra
/// ----------------------------------------------------------------------------
/// rutherford scattering total cross section associated radius approximation
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rs_totra(string const& model) const 
{ 
  double radius(1.);
  if( model == "carli2014" ) { 
    radius = sqrt(rs_totcs("carli2014")/TMath::Pi()); 
  }
  return radius;
}
/// ____________________________________________________________________________





/// dump
/// ----------------------------------------------------------------------------
/// dump properties
void RAtom::dump(ostream &flux) const
{
  flux << "RAtom dump method" << endl; 
  RDump::line(17, "-", flux); 
  flux << "name   : " << getName()                          << endl;
  flux << "symbol : " << getSymbol()                        << endl;
  flux << "Z      : "  << rZ                                << endl;
  flux << "A      : "  << rA                                << endl;
  flux << "mass   : " << getMass()/RUnits::amu << " amu "   << endl;
  flux << "charge : " << getCharge()/RUnits::e   << " e   " << endl;
}
void RAtom::dumpLine(ostream &flux) const
{
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(3); flux << scientific;
  
  flux << setw(8)   << left   << rname    << "   ";
  flux << setw(2)   << right  << rsymbol  << "   ";
  flux << setw(2)   << right  << rZ       << "   ";
  flux << setw(3)   << right  <<  rA      << "   ";
  flux << rmass/RUnits::amu   << "amu"   << "   ";
  flux.precision(1); flux << fixed;
  flux << rcharge/RUnits::e   << "e"     << "   ";
  
  flux << endl; flux.flags(f);  /// restore flags
}

void RAtom::dumpcs(RFlyer* fly, RMachine* mac, std::ostream &flux) const
{
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  int ll(80);
  
  flux << "RAtom dumpcs method ";
  flux << "(rutherford scattering - IPAC2014, C.Carli, TUPRI028)" << endl; 
  RDump::line(ll, "-", flux); 
  flux << "atom:  "; this->dumpLine();
  flux << "flyer: "; fly->dumpLine();
  flux << "       "; fly->dumpBeta();
  flux << "machine: "; mac->dumpLine();
  RDump::line(ll, "-", flux); 
  
  flux.precision(2); flux << scientific;
  flux << "rs_totcs = " << rs_totcs("carli2014")/RUnits::m2 << "m2" << " | ";
  flux << "rs_totra = " << rs_totra()/RUnits::m << "m" << " | ";
  flux << "total cs/radius";
  flux << endl;
  
  flux << "rs_lalcs = " << rs_lalcs(fly,mac)/RUnits::m2 << "m2" << " | ";
  flux << "rs_lalip = " << rs_lalip(fly,mac)/RUnits::m << "m" << " | ";
  flux << "large angle loss cs/impact param";
  flux << endl;
  
  flux << "rs_ebucs = " << rs_ebucs(fly,mac)/RUnits::m2 << "m2" << " | ";
  flux << "rs_ebuip = " << rs_ebuip(fly,mac)/RUnits::m << "m" << " | ";
  flux << "emittance blowup cs/impact param";
  flux << endl;
  
  RDump::line(ll, "-", flux); flux.flags(f);  /// restore flags
}


/// static methods  
/// ****************************************************************************

/// rsip
/// ----------------------------------------------------------------------------
/// rutherford scattering impact parameter for given angle phi
/// IPAC2014, C.Carli, TUPRI028.
double RAtom::rsip(double Z, double z, double m, double v, 
                                   double phi)
{
  double ip(1.);
  /// -> IPAC2014, TUPRI028, Carli
  ip = abs( (Z*z*pow(RUnits::e,2))/(4*M_PI*RUnits::eps0)
            * (2./(m*pow(v,2))) *1./phi 
           );
  /// -> https://en.wikipedia.org/wiki/Rutherford_scattering
  //ip = (z*Z*pow(RUnits::e,2))/(4*M_PI*RUnits::eps0)
  //   * (2./(m*pow(v,2))) *cotan(lossAngle/2.);
  return ip;
}


