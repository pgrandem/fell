/// RFlyer.cc
/// ----------------------------------------------------------------------------
/// rep RFlyer c++ class
/// Pierre Grandemange
/// 19/06/2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <cmath>
/// root classes
/// rep namespaces
#include "RDump.h"
/// rep classes
#include "RFlyer.h"
#include "RParticle.h"
using namespace std;


/// constructors, destructor, copy 
/// ****************************************************************************
RFlyer::RFlyer() : RParticle(), 
  rpos6D{0., 0., 0., 0., 0., 0.}, 
  rtit6D{"x", "y", "z", "px", "py", "pz"}
{}

RFlyer::RFlyer( string const& name, string const& symbol, 
                double charge, double mass, double beta   ) : 
  RParticle(name, symbol, charge, mass), 
  rpos6D{0., 0., 0., 0., 0., betaToMomentum(beta, mass)},
  rtit6D{"x", "y", "z", "px", "py", "pz"}
{}

RFlyer::RFlyer(RParticle* const& par, double beta) : 
  RParticle( par->getName(),par->getSymbol(),par->getCharge(),par->getMass() ),
  rpos6D{0., 0., 0., 0., 0., betaToMomentum(beta, par->getMass())},
  rtit6D{"x", "y", "z", "px", "py", "pz"}
{}

RFlyer::RFlyer(RParticle* const& par, double const x[6]) : 
  RParticle( par->getName(),par->getSymbol(),par->getCharge(),par->getMass() ),
  rpos6D{ x[0], x[1], x[2], x[3], x[4], x[5] },
  rtit6D{"x", "y", "z", "px", "py", "pz"}
{}



RFlyer::~RFlyer()
{}



/// method
/// ****************************************************************************


///particle
/// ----------------------------------------------------------------------------
/// build particle from flyer
RParticle* RFlyer::particle() const 
{ RParticle* par = new RParticle(rname, rsymbol, rcharge, rmass); return par; }
/// set the particle properties of a flyer

void RFlyer::particle(RParticle* par) 
{ 
  rname   = par->getName(); 
  rsymbol = par->getSymbol(); 
  rcharge = par->getCharge(); 
  rmass   = par->getMass();
}
  

/// dump
/// --------------------------------------------------------------------------
/// dump properties
void RFlyer::dump(ostream &flux) const
{
  flux << "RFlyer dump method" << endl; 
  RDump::line(18, "-"); 
  flux << "symbol : " << this->getSymbol() << endl;
  flux << "name   : " << this->getName() << endl;
  flux << "mass   : " << this->getMass()/RUnits::amu      << " amu"     << endl;
  flux << "e0     : " << this->e0()/RUnits::MeV           << " MeV"     << endl;
  flux << "charge : " << this->getCharge()/RUnits::e      << " e"       << endl;
  RDump::line(18, "-"); 
  flux << "beta     : " << this->beta()                                 << endl;
  flux << "gamma    : " << this->gamma()                   							<< endl;
  flux << "ekin     : " << this->ekin()/RUnits::keV        << " kev"    << endl;
  flux << "momentum : " << this->momentum()/RUnits::MeV_c  << " MeV/c"  << endl;
  flux << "velocity : " << this->velocity()/RUnits::m_s    << " m/s"    << endl;
}

void RFlyer::dumpBeta(ostream &flux) const
{
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  
  flux.precision(3); flux << fixed;
  flux << "beta=" 	<< this->beta() 	<< " ";
  flux << "gamma=" 	<< this->gamma() 	<< " ";
  
  flux.precision(3); flux << scientific;
  flux << "ekin=" 	<< this->ekin()/RUnits::MeV << "MeV ";
  flux << "mv=" << this->momentum()/RUnits::MeV_c << "MeV/c ";
  flux << "v=" << this->velocity()/RUnits::m_s    << "m/s ";
  
  flux << endl; flux.flags(f);  /// restore flags
}

void RFlyer::dump6D(ostream &flux) const
{
  int prec(1);        /// precision
  string sep("  ");   /// separator
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(prec); flux << scientific;
  
  flux << "RFlyer dump6D method" << endl; 
  RDump::line(20, "-"); 
  //flux.left(); flux.setw(8);
  flux << right << setw(prec+7) << this->x()/RUnits::mm << "mm" << sep;
  flux << right << setw(prec+7) << this->y()/RUnits::mm << "mm" << sep;
  flux << right << setw(prec+7) << this->z()/RUnits::mm  << "mm" << sep; 
  flux << right << setw(prec+7) << this->px()/RUnits::MeV_c << "MeV_c" << sep; 
  flux << right << setw(prec+7) << this->py()/RUnits::MeV_c << "MeV_c" << sep; 
  flux << right << setw(prec+7) << this->pz()/RUnits::MeV_c << "MeV_c" << sep;
  
  flux << endl; flux.flags(f);  /// restore flags
}

void RFlyer::line6D(ostream &flux) const
{
  /// save previous stream format settings
  ios::fmtflags f(flux.flags());        /// save current flux flags
  streamsize    p(flux.precision());    /// save current flux precision
  /// set new stream format settings
  flux.flags(ios::right | ios::scientific); /// ajust on right, scientific
  int prec(1);                              /// precision
  flux.precision(prec);                     /// set precision 
  int fs(7);            /// format size
  int w(prec+fs);       /// width
  string sep(" ");     /// separator
  
  /// debug
  //double tst(31415.16171819);
  //flux << f << sep;
  //flux << p << sep;
  
  /// stream output
  flux << setw(w) << x()/RUnits::mm << RUnits::mm.getSymbol() << sep;
  flux << setw(w) << y()/RUnits::mm << RUnits::mm.getSymbol() << sep;
  flux << setw(w) << z()/RUnits::mm << RUnits::mm.getSymbol() << sep; 
  flux << setw(w) << px()/RUnits::MeV_c << RUnits::MeV_c.getSymbol() << sep; 
  flux << setw(w) << py()/RUnits::MeV_c << RUnits::MeV_c.getSymbol() << sep; 
  flux << setw(w) << pz()/RUnits::MeV_c << RUnits::MeV_c.getSymbol() << sep;
  flux << endl;
  
  /// debug 
  /*flux << tst << sep;
  flux.precision(4);
  flux << tst << sep;
  flux.precision(10);
  flux << tst << sep;
  flux.unsetf(ios_base::floatfield);
  flux << tst << sep;
  flux.flags(f);        /// restore flags
  flux.precision(p);
  flux << tst << sep;
  flux << endl; */
  
  /// restore previous stream format settings
  flux.flags(f);        /// restore flags
  flux.precision(p);    /// restore precision
  
}



/// static methods
/// ****************************************************************************

/// betaToEkin
/// ----------------------------------------------------------------------------
/// ekin = [gamma(beta)-1]*e0
double RFlyer::betaToEkin(double beta, double mass)
{ return (RFlyer::betaToGamma(beta)-1)*massToE0(mass); }

/// betaToGamma
/// ----------------------------------------------------------------------------
/// gamma = 1./sqrt(1-beta2)
double RFlyer::betaToGamma(double beta) 
{ return 1./sqrt(1.-pow(beta, 2)); }

/// betaToMomentum
/// ----------------------------------------------------------------------------
/// momentum = gamma(beta)*m0*beta*c
double RFlyer::betaToMomentum(double beta, double m0) 
{ return (betaToGamma(beta)*m0*beta*RUnits::c); }

/// betaToVelocity
/// ----------------------------------------------------------------------------
/// velocity = beta*c
double RFlyer::betaToVelocity(double beta) 
{ return beta*RUnits::c; }

/// ekinToBeta
/// ----------------------------------------------------------------------------
/// beta = gammaToBeta(ekinToGamma)
double RFlyer::ekinToBeta(double ekin, double mass) 
{ return gammaToBeta(ekinToGamma(ekin, mass)); }

/// ekinToGamma
/// ----------------------------------------------------------------------------
/// gamma = ekin/e0(mass)+1
double RFlyer::ekinToGamma(double ekin, double mass) 
{ return ekin/massToE0(mass)+1; }

/// ekinToMomentum
/// ----------------------------------------------------------------------------
/// momentum = gamma*beta*m0*c
double RFlyer::ekinToMomentum(double ek, double m0) 
{ return ekinToGamma(ek, m0) * ekinToBeta(ek, m0) * m0 * RUnits::c; }

/// gammaToBeta
/// ----------------------------------------------------------------------------
/// beta = sqrt(1-1/gamma2)
double RFlyer::gammaToBeta(double gamma)
{ return sqrt(1.-1./pow(gamma, 2)); }

/// gammaToEkin
/// ----------------------------------------------------------------------------
/// ekin = (gamma-1)*e0
double RFlyer::gammaToEkin(double gamma, double mass)
{ return (gamma-1.)*massToE0(mass); }

/// massToE0
/// ----------------------------------------------------------------------------
/// e0 = m*c2 !!! ;)
double RFlyer::massToE0(double mass) 
{ return mass*RUnits::c*RUnits::c; }

/// momentumToBeta
/// ----------------------------------------------------------------------------
/// beta = p * 1./sqrt( (m0*c)2 + p2 )
double RFlyer::momentumToBeta(double p, double m0) 
{ return p/sqrt( pow(m0*RUnits::c,2) + pow(p,2) ); }

/// momentumToEkin
/// ----------------------------------------------------------------------------
/// Ekin = ( gamma(p)-1 ) * E0
double RFlyer::momentumToEkin(double p, double m0) 
{ return (betaToGamma(momentumToBeta(p, m0))-1)* massToE0(m0); }

