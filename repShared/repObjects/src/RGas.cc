/// RGas.cc
/// ----------------------------------------------------------------------------
/// rep RGas c++ class 
/// Pierre GRANDEMANGE
/// 15/06/2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
/// root classes
#include "TMath.h"
/// rep classes
#include "RBeam.h"
#include "RFlyer.h"
#include "RMachine.h"
#include "RMolecule.h"
#include "RObject.h"
#include "RGas.h"
/// rep namespace
#include "RDump.h"
#include "RUnits.h"

using namespace std;

/// constructor, destructor, copy 
/// **************************************************************************
RGas::RGas(std::string const& name) : RObject(name), 
  rcompound(), rfraction(), rpressureN2eq(0), rtemperature(0),
  rgauge("pkr261")
{}

RGas::RGas(std::string const& name, double pN2eq, double temperature) : 
  RObject(name), rcompound(), rfraction(), rpressureN2eq(pN2eq), 
  rtemperature(temperature), rgauge("pkr261")
{}

RGas::RGas(RMolecule* const& mol, double pN2eq, double T) : 
  RObject(mol->getName()), rfraction(1, 1.), rpressureN2eq(pN2eq), 
  rtemperature(T), rgauge("pkr261")
{
  rcompound.push_back(mol);
}


RGas::~RGas()
{
  /// delete molecule pointers in compound vector
  for( int i=0; i<rcompound.size(); ++i) { delete rcompound[i]; rcompound[i]=0;}
}





/// dump methods
/// **************************************************************************

/// dump
/// ----------------------------------------------------------------------------
/// dump properties
void RGas::dump(ostream& flux) const
{
  ios::fmtflags f(flux.flags());    /// save current flux flags
  streamsize p(flux.precision());   /// save current flux precision
  
  double ll(67);  /// line length;
  flux << "RGas dump method" << endl; 
  RDump::line(16, "-", flux); 
  vector<double> den(this->density(true)); /// show pK factor warning
  
  /// start array
  flux << left << setw(15) << rname << " | ";
  
  flux.precision(1); flux << scientific;
  flux << "p=" << rpressureN2eq/RUnits::mbar << "mbar_N2eq" << "  ";
  
  flux.precision(1); flux << fixed;
  flux << "T=" << rtemperature/RUnits::K << "K" << "  ";
  
  flux.precision(3); flux << scientific;
  flux << "n=" << (this->densityN2eq())*RUnits::m3 << "/m3_N2eq" << "  ";
  
  flux << endl; flux.flags(f);  /// restore flags
  //RDump::line(ll, "-", flux); flux.flags(f);  /// restore flags
  
  for( int j=0; j<rcompound.size(); ++j) { 
    ios::fmtflags f(flux.flags());    /// save current flags in flux
  	
  	RDump::line(ll, "-", flux);
    flux << (rfraction[j])/RUnits::percent << "%" << " "; 
    flux << (rcompound[j]->getSymbol()) << "("; 
    flux << (rcompound[j]->getSymbolIsotope()) << ") | " ;
    
    flux.precision(2); flux << scientific;
    flux << "n=" << den[j]*RUnits::m3 << "/m3(pK=";
    flux.precision(1); flux << fixed;
    flux << (rcompound[j]->pKfactor(rpressureN2eq, rgauge, false)) << ")";
    flux << " | ";
    
    //flux << (rcompound[j]->getName()) << " " ;
    //flux << (rcompound[j]->getSymbolIsotope()) << " " ;
    flux.precision(2); flux << fixed;
    flux << "m=" << (rcompound[j]->getMass())/RUnits::amu << "amu" << "   " ;
    //flux << (rcompound[j]->getCharge())/RUnits::e << " e" << "   " ;
    //flux << (rcompound[j]->rs_totcs()) << " m^2" << "   ";
    
    flux << endl; flux.flags(f);  /// restore flags
    //RDump::line(ll, "-", flux); flux.flags(f);  /// restore flags
  
    for( int k=0; k<(rcompound[j]->getComposition()).size(); ++k) { 
      ((rcompound[j]->getComposition())[k])->dumpLine(flux); 
    }
  }
  RDump::line(ll, "-", flux); 
  flux.flags(f);      /// restore flags
  flux.precision(p);  /// restore precision
}



/// dumpLine
/// ----------------------------------------------------------------------------
void RGas::dumpLine(ostream &flux) const
{
	ios::fmtflags f(flux.flags());	/// save current flags in flux
  
  flux << this->getName() << " | ";
  flux.precision(2); flux << scientific;
  flux << "pN2eq=" << this->getPressureN2eq()/RUnits::mbar << "mbar - ";;
  flux.precision(1); flux << fixed;
  flux << "T=" << this->getTemperature()/RUnits::K << "K - ";
  flux.precision(2); flux << scientific;
  flux << "nN2eq=" << this->densityN2eq()/RUnits::m3 << "/m3";;
  
  flux << endl; flux.flags(f);		/// restore flags
}


/// dumpRates (flyer version)
/// ----------------------------------------------------------------------------
void RGas::dumpRates(RFlyer* const& fly, RMachine* const & mac, 
                     ostream& flux) const
{
  /// header
  int ll(80); /// line length
  string sep("  "); /// default separator
  flux << "RGas dumpRates method (with flyer)" << endl; RDump::line(34, "-");
  
  /// usefull parameters
  int size(this->getCompound().size());  ///nb of molecule
  vector<double> ics(this->ion_pbarcs(fly));  /// cross section vectors
  vector<double> tcs(this->rs_totcs());
  vector<double> lcs(this->rs_lalcs(fly, mac));
  vector<double> ecs(this->rs_ebucs(fly, mac));
  vector<double> ira(this->ion_pbarnu(fly));    /// interaction rate vectors
  vector<double> tra(this->rs_totnu(fly));
  vector<double> lra(this->rs_lalnu(fly, mac));
  vector<double> era(this->rs_ebunu(fly, mac));
  vector<double> den(this->density(true));  /// vec of densities
   
  /// gas - flyer - machine
  cout << "gas: "; 			this->dumpLine(flux);
  cout << "flyer: "; 		fly->dumpLine(flux);
  cout << "       "; 		fly->dumpBeta(flux);
  cout << "machine: "; 	mac->dumpLine(flux);
  
  /// dump gas rates
  RDump::line(ll, "-", flux);
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(3); flux << scientific; flux << left;
  /// first line
  flux << setw(4) 	<< " " << sep;
  flux << setw(14) 	<< "ionisation" << sep;
  flux << setw(14) 	<< "scatt tot" << sep;
  flux << setw(14) 	<< "las loss" << sep;
  flux << setw(14) 	<< "em blow up" << sep;
  flux << endl;
  flux << setw(4) << "rate" << sep;
  flux << setw(9) << ion_pbarnuGas(fly) 		<< setw(5) << "/s" << sep;
  flux << setw(9) << rs_totnuGas(fly) 			<< setw(5) << "/s" << sep;
  flux << setw(9) << rs_lalnuGas(fly, mac) 	<< setw(5) << "/s" << sep;
  flux << setw(9) << rs_ebunuGas(fly, mac)	<< setw(5) << "/s" << sep;
  flux << endl;
  
  flux << fixed;
  flux.flags(f);  /// restore flags
  
  /// detail molecules rates and cs
  RDump::line(ll, "-", flux);
  for( int i=0; i<this->getCompound().size(); ++i ) {
    flux << (rfraction[i])/RUnits::percent << "%" << " "; 
    flux << (rcompound[i]->getSymbol()) << "  "; 
    flux << "n=" << den[i]*RUnits::m3 << "/m3(pK=";
    flux << (rcompound[i]->pKfactor(rpressureN2eq, rgauge, false)) << ")";
    flux << " | ";
    flux << endl;
    
    flux.precision(3); 
    flux << scientific;
  
    ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(3); flux << scientific; flux << left;
  /// first line
  flux << setw(4) 	<< " " << sep;
  flux << setw(14) 	<< "ionisation" << sep;
  flux << setw(14) 	<< "scatt tot" << sep;
  flux << setw(14) 	<< "las loss" << sep;
  flux << setw(14) 	<< "em blow up" << sep;
  flux << endl;
  flux << setw(4) << "cs" << sep;
  flux << setw(9) << ics[i] << setw(5) << "m2" << sep;
  flux << setw(9) << tcs[i]	<< setw(5) << "m2" << sep;
  flux << setw(9) << lcs[i]	<< setw(5) << "m2" << sep;
  flux << setw(9) << ecs[i]	<< setw(5) << "m2" << sep;
  flux << endl;
  flux << setw(4) << "rate" << sep;
  flux << setw(9) << ira[i] << setw(5) << "/s" << sep;
  flux << setw(9) << tra[i]	<< setw(5) << "/s" << sep;
  flux << setw(9) << lra[i]	<< setw(5) << "/s" << sep;
  flux << setw(9) << era[i]	<< setw(5) << "/s" << sep;
  flux << endl;
  
  flux << fixed;
  flux.flags(f);  /// restore flags
  
  }
  
  /// end of array line
  RDump::line(ll, "-", flux);
  flux.flags(f);  /// restore flags
}



/// dumpRates (beam version)
/// ----------------------------------------------------------------------------
void RGas::dumpRates(RBeam* const& bea, RMachine* const & mac, 
                     ostream& flux) const
{
  /// header
  int ll(80); /// line length
  string sep("  "); /// default separator
  flux << "RGas dumpRates method (for beam)" << endl; RDump::line(32, "-");
  
  /// usefull parameters
  RFlyer*        fly( bea->getFlyer()               );
  int            size( this->getCompound().size()   ); /// nb of molecule
  vector<double> ics( this->ion_pbarcs(fly)         ); /// cross section vectors
  vector<double> tcs( this->rs_totcs()              );
  vector<double> lcs( this->rs_lalcs(fly, mac)      );
  vector<double> ecs( this->rs_ebucs(fly, mac)      );
  vector<double> ira( this->ion_pbarnu(fly)         ); /// inter rate vectors
  vector<double> tra( this->rs_totnu(fly)           );
  vector<double> lra( this->rs_lalnu(fly, mac)      );
  vector<double> era( this->rs_ebunu(fly, mac)      );
  vector<double> den( this->density(true)           );  /// vec of densities
   
  /// gas - flyer - machine
  cout << "gas: "; 			this->dumpLine(flux);
  cout << "flyer: "; 		fly->dumpLine(flux);
  cout << "       "; 		fly->dumpBeta(flux);
  cout << "machine: "; 	mac->dumpLine(flux);
  
  /// dump gas rates
  RDump::line(ll, "-", flux);
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(3); flux << scientific; flux << left;
  /// first line
  flux << setw(14) 	<< "rate unit"  << sep;
  flux << setw(14) 	<< "ionisation" << sep;
  flux << setw(14) 	<< "scatt tot"  << sep;
  flux << setw(14) 	<< "las loss"   << sep;
  flux << setw(14) 	<< "em blow up" << sep;
  flux << endl;
  flux << setw(14) << "/ion/s" << sep;
  flux << setw(14) << ion_pbarnuGas(fly) 		<< sep;
  flux << setw(14) << rs_totnuGas(fly) 			<< sep;
  flux << setw(14) << rs_lalnuGas(fly, mac) << sep;
  flux << setw(14) << rs_ebunuGas(fly, mac) << sep;
  flux << endl;
  flux << setw(14) << "/beam/m/s" << sep;
  flux << setw(14) << ion_pbarnuGas(bea, mac)*RUnits::m*RUnits::s << sep;
  flux << setw(14) << rs_totnuGas(bea, mac)*RUnits::m*RUnits::s << sep;
  flux << setw(14) << rs_lalnuGas(bea, mac)*RUnits::m*RUnits::s << sep;
  flux << setw(14) << rs_ebunuGas(bea, mac)*RUnits::m*RUnits::s << sep;
  flux << endl;
  flux << setw(14) << "/beam/cm/100ms" << sep;
  flux << setw(14) << ion_pbarnuGas(bea, mac)*RUnits::cm*100*RUnits::ms << sep;;
  flux << setw(14) << rs_totnuGas(bea, mac)*RUnits::cm*100*RUnits::ms << sep;; 
  flux << setw(14) << rs_lalnuGas(bea, mac)*RUnits::cm*100*RUnits::ms << sep;; 
  flux << setw(14) << rs_ebunuGas(bea, mac)*RUnits::cm*100*RUnits::ms << sep;; 
  flux << endl;
  
  flux << fixed;
  flux.flags(f);  /// restore flags
  
  /// detail molecules rates and cs
  RDump::line(ll, "-", flux);
  for( int i=0; i<this->getCompound().size(); ++i ) {
    flux << (rfraction[i])/RUnits::percent << "%" << " "; 
    flux << (rcompound[i]->getSymbol()) << "  "; 
    flux << "n=" << den[i]*RUnits::m3 << "/m3(pK=";
    flux << (rcompound[i]->pKfactor(rpressureN2eq, rgauge, false)) << ")";
    flux << " | ";
    flux << endl;
    
    flux.precision(3); 
    flux << scientific;
  
    ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(3); flux << scientific; flux << left;
  /// first line
  flux << setw(14) 	<< "val|unit" << sep;
  flux << setw(14) 	<< "ionisation" << sep;
  flux << setw(14) 	<< "scatt tot" << sep;
  flux << setw(14) 	<< "las loss" << sep;
  flux << setw(14) 	<< "em blow up" << sep;
  flux << endl;
  flux << setw(14) << "cs|m2" << sep;
  flux << setw(14) << ics[i] 	<< sep;
  flux << setw(14) << tcs[i]	<< sep;
  flux << setw(14) << lcs[i]	<< sep;
  flux << setw(14) << ecs[i]	<< sep;
  flux << endl;
  flux << setw(14) << "rate|/ion/s" << sep;
  flux << setw(14) << ira[i] 	<< sep;
  flux << setw(14) << tra[i]	<< sep;
  flux << setw(14) << lra[i]	<< sep;
  flux << setw(14) << era[i]	<< sep;
  flux << endl;
  
  flux << fixed;
  flux.flags(f);  /// restore flags
  
  }
  
  /// end of array line
  RDump::line(ll, "-", flux);
  flux.flags(f);  /// restore flags
}

















/// "other" methods
/// **************************************************************************

/// checkComp
/// --------------------------------------------------------------------------
/// check if sum of fractional molecular componants = 1;
bool RGas::checkComp() const
{
  bool check=true;
  double sum(0.);
  for(int i=0; i<rfraction.size(); ++i) { sum += rfraction[i]; }
  if ( sum!=1. || rcompound.size()!=rfraction.size() ) {
    cout << endl; 
    cout << "  rep WARNING : gas molecule size and/or fraction parts" << endl;
    cout << "  rfracSum = " << sum << endl;
    cout << "  rfraction.size() = " << rfraction.size() << endl;
    cout << "  rcompound.size() = " << rfraction.size() << endl;
    check=false;
    
  }
  return check;
}
  
/// density
/// --------------------------------------------------------------------------
/// Compute density for each molecule in the gas taking fractionnal
/// composition and pK correction factor into account.
vector<double> RGas::density(bool verbose) const
{
  vector<double> ppv(this->pPartial(verbose));   /// call checkComp
  vector<double> vec(rcompound.size(), 0);
  for( int i=0; i<vec.size(); ++i ) { 
    vec[i] = ppv[i]/(RUnits::kB*rtemperature);
  }
  return vec;
}

/// densityFracN2eq
/// ----------------------------------------------------------------------------
/// Compute N2eq density for each molecule in the gas taking fractionnal
/// composition into account.
vector<double> RGas::densityFracN2eq() const 
{
  vector<double> vec(0);
  if( checkComp() ) { /// no fractionnal problem
    for( int i=0; i<rcompound.size(); ++i ) {
      vec.push_back(rfraction[i]*densityN2eq()); /// pi = fractional * pN2eq
    }
  }
  return vec;
}


/// pPartialN2eq
/// ----------------------------------------------------------------------------
/// get N2eq partial pressures of molecule constituating the gas
/// Not physical
vector<double> RGas::pPartialN2eq() const
{
  vector<double> piv(0);
  if( checkComp() ) { /// no fractionnal problem
    for( int i=0; i<rcompound.size(); ++i ) {
      piv.push_back(rfraction[i]*rpressureN2eq); /// pi = fractional * pN2eq
    }
  }
  return piv;
}

/// pPartial
/// --------------------------------------------------------------------------
/// get partial pressures of molecule constituating the gas
/// taking into account fractional and correction factor for molecules
vector<double> RGas::pPartial(bool verbose) const
{
  vector<double> pn2(this->pPartialN2eq()); /// get pN2eq vec, checkComp
  vector<double> vec(0);
  for( int i=0; i<rcompound.size(); ++i ) {
    vec.push_back(
      (rcompound[i]->pKfactor(rpressureN2eq,rgauge,verbose))*pn2[i] );
  }
  return vec;
}




/// cross sections, interaction rates and gas rates
/// ----------------------------------------------------------------------------
  
/// cross sections per molecules
/// ----------------------------------------------------------------------------
vector<double> RGas::ion_pbarcs(RFlyer* const& fly) const
{
  vector<double> vec;
  for( int i=0; i<this->getCompound().size(); ++i ) {
    vec.push_back(this->getCompound()[i]->ion_pbarcs(fly));
  }
  return vec;
}

vector<double> RGas::rs_ebucs(RFlyer* const& fly, RMachine* const& mac) const
{
  vector<double> vec;
  for( int i=0; i<this->getCompound().size(); ++i ) {
    vec.push_back(this->getCompound()[i]->rs_ebucs(fly, mac));
  }
  return vec;
}

vector<double> RGas::rs_lalcs(RFlyer* const& fly, RMachine* const& mac) const
{
  vector<double> vec;
  for( int i=0; i<this->getCompound().size(); ++i ) {
    vec.push_back(this->getCompound()[i]->rs_lalcs(fly, mac));
  }
  return vec;
}

vector<double> RGas::rs_totcs() const
{
  vector<double> vec;
  for( int i=0; i<this->getCompound().size(); ++i ) {
    vec.push_back(this->getCompound()[i]->rs_totcs());
  }
  return vec;
}


/// interaction rates per molecules (per flyer)
/// ----------------------------------------------------------------------------
vector<double> RGas::ion_pbarnu(RFlyer* const& fly) const
{
  vector<double> vec(this->getCompound().size()); /// output
  vector<double> ics = ion_pbarcs(fly);           /// interaction cross section
  vector<double> den = density(false);                 /// molecule density
  
  for( int i=0; i<vec.size(); ++i ) { 
    vec[i] = ics[i] * den[i] * (fly->beta()) * RUnits::c;
    //cout << ics[i] << "   " << den[i] << "   " << (fly->getBeta()) << "   ";
    //cout << RUnits::c << "   " << vec[i] << endl;
    
  }
  return vec;
}

vector<double> RGas::rs_ebunu(RFlyer* const& fly, RMachine* const& mac) const
{
  vector<double> vec(this->getCompound().size()); /// output
  vector<double> ics = rs_ebucs(fly, mac);        /// interaction cross section
  vector<double> den = density(false);         /// molecule density
  
  for( int i=0; i<vec.size(); ++i ) { 
    vec[i] = ics[i] * den[i] * (fly->beta()) * RUnits::c;
    //cout << ics[i] << "   " << den[i] << "   " << (fly->getBeta()) << "   ";
    //cout << RUnits::c << "   " << vec[i] << endl;
  }
  return vec;
}

vector<double> RGas::rs_ebunuCarli(RFlyer* const& fly, 
                                   RMachine* const& mac) const
{
  vector<double> vec(rs_ebunu(fly, mac));    /// output
  for( int i=0; i<vec.size(); ++i ) { 
    vec[i] *= mac->getAcceptance();
  }
  return vec;
}

vector<double> RGas::rs_lalnu(RFlyer* const& fly, RMachine* const& mac) const
{
  vector<double> vec(this->getCompound().size()); /// output
  vector<double> ics = rs_lalcs(fly, mac);        /// interaction cross section
  vector<double> den = density(false);         /// molecule density
  
  for( int i=0; i<vec.size(); ++i ) { 
    vec[i] = ics[i] * den[i] * (fly->beta()) * RUnits::c;
    //cout << ics[i] << "   " << den[i] << "   " << (fly->getBeta()) << "   ";
    //cout << RUnits::c << "   " << vec[i] << endl;
  }
  return vec;
}

vector<double> RGas::rs_totnu(RFlyer* const& fly) const
{
  vector<double> vec(this->getCompound().size()); /// output
  vector<double> ics = rs_totcs();        /// interaction cross section
  vector<double> den = density(false);         /// molecule density
  
  for( int i=0; i<vec.size(); ++i ) { 
    vec[i] = ics[i] * den[i] * (fly->beta()) * RUnits::c;
    //cout << ics[i] << "   " << den[i] << "   " << (fly->getBeta()) << "   ";
    //cout << RUnits::c << "   " << vec[i] << endl;
  }
  return vec;
}



/// rates for the gas per beam (times number of flyers, divided by mac length)
/// ----------------------------------------------------------------------------
double RGas::ion_pbarnuGas(RFlyer* const& fly) const
{
  double dbl(0); vector<double> vec(ion_pbarnu(fly));
  for( int i=0; i<vec.size(); ++i ) { dbl += vec[i]; }
  return dbl;
}

double RGas::rs_ebunuGas(RFlyer* const& fly,RMachine* const& mac) const
{
  double dbl(0); vector<double> vec(rs_ebunu(fly, mac));
  for( int i=0; i<vec.size(); ++i ) { dbl += vec[i]; }
  return dbl;
}

double RGas::rs_lalnuGas(RFlyer* const& fly,RMachine* const& mac) const
{
  double dbl(0); vector<double> vec(rs_lalnu(fly, mac));
  for( int i=0; i<vec.size(); ++i ) { dbl += vec[i]; }
  return dbl;
}

double RGas::rs_totnuGas(RFlyer* const& fly) const
{
  double dbl(0); vector<double> vec(rs_totnu(fly));
  for( int i=0; i<vec.size(); ++i ) { dbl += vec[i]; }
  return dbl;
}
  

/// rates for the gas per beam (times number of flyers, divided by mac length)
/// 2017 09 20 : to be change to adapt to new flyer/beam classes
/// ----------------------------------------------------------------------------
double RGas::ion_pbarnuGas(RBeam* const& beam,  RMachine* const& mac) const
{
	RFlyer* fly( beam->getFlyer() );
	double 	npb( beam->getN()     );
	double 	len( mac->getLength() );
	
	return ion_pbarnuGas(fly) * npb/len;
}

double RGas::rs_ebunuGas(RBeam* const& beam, RMachine* const& mac) const
{
	RFlyer* fly( beam->getFlyer() );
	double 	npb( beam->getN()     );
	double 	len( mac->getLength() );
	
	return rs_ebunuGas(fly, mac) * npb/len;
}

double RGas::rs_lalnuGas(RBeam* const& beam, RMachine* const& mac) const
{
	RFlyer* fly( beam->getFlyer() );
	double 	npb( beam->getN()     );
	double 	len( mac->getLength() );
	
	return rs_lalnuGas(fly, mac) * npb/len;
}

double RGas::rs_totnuGas(RBeam* const& beam,  RMachine* const& mac) const
{
	RFlyer* fly( beam->getFlyer() );
	double 	npb( beam->getN()     );
	double 	len( mac->getLength() );
	
	return rs_totnuGas(fly) * npb/len;
}



/// rsRates()
/// --------------------------------------------------------------------------
/// computes beam/gas mean interaction rates and standard deviation(sigma)
/// sum flyer/gas interaction rates over the beam flyer collection
/// two computation loop to get these values
/// returns a ptr:
///   - ptr[0] | ion_pbarnuGas(beam) = sum( ion_pbarnuGas(flyer)  )/macLength
///   - ptr[1] | rs_ebunuGas(beam)   = sum( rs_ebunuGas(flyer)    )/macLength
///   - ptr[2] | rs_lalnuGas(beam)   = sum( rs_lalnuGas(flyer)    )/macLength
///   - ptr[3] | rs_totnuGas(beam)   = sum( rs_totnuGas(flyer)    )/macLength
///   - ptr[4] | sigma( ion_pbarnuGas(beam) )
///   - ptr[5] | sigma( rs_ebunuGas(beam)   )
///   - ptr[6] | sigma( rs_lalnuGas(beam)   )
///   - ptr[7] | sigma( rs_totnuGas(beam)   )
double* RGas::rsRates(RBeam* const& beam, RMachine* const& mac) const
{
  ///parameters
  double* rat= new double[8]{0.};             /// output
  RFlyer* fco( beam->getFlyColl() );  /// flyer collection
  int     nfl( beam->getN()  );               /// number of flyers
  double 	mle( mac->getLength() );            /// machine length
  
  /// 1st nfl loop : sum of xi
  for( int i=0; i<beam->getN(); ++i ) { 
    //rat[0] += ion_pbarnuGas( (RFlyer*)&fco[i]      );
    rat[0] += 1.;
    rat[1] += rs_ebunuGas(   (RFlyer*)&fco[i], mac );
    rat[2] += rs_lalnuGas(   (RFlyer*)&fco[i], mac );
    rat[3] += rs_totnuGas(   (RFlyer*)&fco[i]      );
  }
  /// flyer coll nu mean value
  for( int i=0; i<4; ++i ) { rat[i] /= nfl; }
  
  /// 2nd nfl loop : sum of (xi - mean)^2
  for( int i=0; i<beam->getN(); ++i ) { 
    //rat[4] += pow( (ion_pbarnuGas(&fco[i]) - rat[0]),    2 );
    rat[4] += 1.;
    rat[5] += pow( (rs_ebunuGas(&fco[i], mac) - rat[1]), 2 );
    rat[6] += pow( (rs_lalnuGas(&fco[i], mac) - rat[2]), 2 );
    rat[7] += pow( (rs_totnuGas(&fco[i]) - rat[3]),      2 );
  }
  /// flyer coll nu sigma
  for( int i=4; i<8; ++i ) { rat[i] = sqrt(1./nfl*rat[i]); }
  return rat;
}








