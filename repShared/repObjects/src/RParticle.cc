/// RParticle.cc
/// ----------------------------------------------------------------------------
/// rep RParticle c++ class
/// Pierre Grandemange
/// 12/06/2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// rep namespaces
#include "RUnits.h"
#include "RDump.h"
/// rep classes
#include "RParticle.h"
using namespace std;



/// static attributes initialisation
/// ****************************************************************************
int RParticle::rnb(0);                          /// nb of running instances
vector<string> RParticle::rlist(1,"emptyList"); /// list of running instances



/// constructors and destructor
/// ****************************************************************************
RParticle::RParticle() 
  : RObject("parName"), rcharge(0.), rmass(0.), rsymbol("parSymb")
{ 
  ++rnb;
  if( rlist[0] == "emptyList" ) { rlist[0]=rname; }
  else                          { rlist.push_back(rname); }
}

RParticle::RParticle(string name, string symbol) 
  : RObject(name), rcharge(0.), rmass(0.), rsymbol(symbol)
{ 
  ++rnb;
  if( rlist[0] == "emptyList" ) { rlist[0]=rname; }
  else                          { rlist.push_back(rname); }
}

RParticle::RParticle(string name, string symbol, double charge,double mass)
  : RObject(name), rcharge(charge), rmass(mass), rsymbol(symbol)
{ 
  ++rnb;
  if( rlist[0] == "emptyList" ) { rlist[0]=rname; }
  else                          { rlist.push_back(rname); }
}

RParticle::~RParticle() 
{ 
  --rnb;
}



/// methods
/// ****************************************************************************
/// dump properties
/// ----------------------------------------------------------------------------
void RParticle::dump(ostream &flux) const
{
  flux << "RParticle dump method" << endl; 
  RDump::line(21, "-", flux); 
  flux << "name : "    << getName() << endl; 
  flux <<  "symbol : " << getSymbol() << endl;
  flux <<  "mass   : " << getMass()/RUnits::MeV_c2 << " MeV/c2" << endl;
  flux <<  "charge : " << getCharge()/RUnits::e << " e"    << endl;
}
void RParticle::dumpLine(ostream &flux) const
{
	ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(3); flux << fixed; flux << left;
  
  flux << rname             					<< ": ";
  flux << rsymbol           					<< " ";
  flux << "m=" << rmass/RUnits::amu << "amu"  << " ";
  
  flux.precision(1);
  flux << "C=" << rcharge/RUnits::e << "e"    << " ";
  
  flux << endl; flux.flags(f);  /// restore flags
}


/// static methods
/// ****************************************************************************
int RParticle::nb() { return rnb;}

void RParticle::list(ostream &flux)
{
  flux << "RParticle list method" << endl; 
  RDump::line(21, "-", flux); 
  cout << "Number of running instances :  " << rnb << endl;
  if( rnb==0 ) {
    cout << "the list is empty" << endl;
  } 
  else {
    cout << "list of particles instantiated (even if deleted) : " ;
    cout << rlist.size() << endl;
    for( int i=0; i<rlist.size(); ++i ) {
      cout << setw(2) << i << "  " << rlist[i] << endl;
    }
  }
}



/// external operators
/// ****************************************************************************
/*ostream& operator<<(ostream &flux, RParticle const& particle)
{
  particle.dump(flux);
  return flux;
}*/



