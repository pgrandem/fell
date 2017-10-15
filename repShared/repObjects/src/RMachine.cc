/// RMachine.h
/// ----------------------------------------------------------------------------
/// rep RMachine c++ class 
/// Pierre GRANDEMANGE
/// 19/06/2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
/// root classes
//#include "TCanvas.h"
/// rep namespaces
#include "RUnits.h"
/// rep classes
#include "RMachine.h"
#include "RObject.h"
using namespace std;


/// constructors, destructor, copy 
/// ****************************************************************************
RMachine::RMachine() 
  : RObject(), racceptance(1.), rlength(1.)
{
	rtwissBeta[0]=1.;
	rtwissBeta[1]=2.;
}

RMachine::RMachine(string const& name) 
  : RObject(name), racceptance(1.), rlength(1.)
{
	rtwissBeta[0]=1.;
	rtwissBeta[1]=2.;
}


RMachine::RMachine(string const& name, double length) 
  : RObject(name), racceptance(1.), rlength(length)
{
	rtwissBeta[0]=1.;
	rtwissBeta[1]=2.;
}

RMachine::RMachine(string const& name, double length, 
									 double accept, double betaT) : RObject(name), 
	racceptance(accept), rlength(length)
{
	rtwissBeta[0]=betaT;
	rtwissBeta[1]=betaT;
}

RMachine::~RMachine()
{}



/// methods
/// ****************************************************************************

/// dump
/// ----------------------------------------------------------------------------
/// dump properties
void RMachine::dumpLine(ostream &flux) const
{
  ios::fmtflags f(flux.flags());    /// save current flags in flux
  flux.precision(6); ///flux << fixed;
  
  flux << this->getName() << ": ";
  flux << "L=" 			  << this->getLength()/RUnits::m << "m | ";
  flux << "tBetaX=" 	<< this->getTwissBetaX()/RUnits::m << "m | ";
  flux << "tBetaY=" 	<< this->getTwissBetaY()/RUnits::m << "m | ";
  flux << "A=" 			  << this->getAcceptance()/RUnits::um << "um | ";
  flux << "Axangle="  << this->acceptanceXangle() << " | ";
  flux << "Ayangle="  << this->acceptanceYangle() << "";
  
  flux << endl; flux.flags(f);  /// restore flags
}
