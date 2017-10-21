/// RDump.cc
/// ----------------------------------------------------------------------------
/// rep RDump c++ namespace
/// 23/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <string>
/// rep classes
/// rep namespaces
#include "RDump.h"
using namespace std;


namespace RDump {

/// line
/// --------------------------------------------------------------------------
/// dump a line
void line(int width, string const& charac, ostream& flux)
{ for( int i=0; i<width; ++i ) { flux << charac; } flux << endl; }


} /// end of namespace
