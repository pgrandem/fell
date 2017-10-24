/// RObject.cc
/// ----------------------------------------------------------------------------
/// rep RObject c++ object
/// original: scrapers/TestEm5 geant4 analysis program
/// 12/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <string>
/// root librairies
/// rep objects
#include "RObject.h"
/// rep namespaces
#include "RDump.h"
/// rep shared functions
using namespace std;



/// static attributes initialisation
/// ****************************************************************************
int RObject::rnb(0);                            /// nb of running instances


// constructor, destructor, copy
/// ****************************************************************************
RObject::RObject() : rname("objName") { ++rnb; }
RObject::RObject(const string& name) : rname(name) { ++rnb; }
RObject::~RObject() { --rnb; }
RObject::RObject(RObject const& objToCopy)
{
  rname = objToCopy.rname;
  { ++rnb; }
}



/// methods
/// ****************************************************************************

/// dump properties
/// ----------------------------------------------------------------------------
void RObject::dump(ostream &flux) const
{
  RDump::header("RObject::dump \"" + this->getName() + "\"");
  //RDump::line(19, "-", flux); 
  flux << "name : " << this->getName() << endl; 
}



/// operators
/// ****************************************************************************

/// = operator
/// ----------------------------------------------------------------------------
RObject& RObject::operator=(const RObject& otherObj)
{
  if( this != &otherObj ) {
    rname = otherObj.rname;
  }
  return *this;
}

/// < operator
/// ----------------------------------------------------------------------------
bool RObject::operator<(const RObject& otherObj) const
{
  //cout << "you're in!" << endl;
  //return 0;
  return rname < otherObj.rname;
}



/// static methods
/// ****************************************************************************
int RObject::nb() { return rnb;}



/// external operators
/// ****************************************************************************
ostream& operator<<(ostream &flux, RObject const& object)
{
  object.dump(flux);
  return flux;
}





