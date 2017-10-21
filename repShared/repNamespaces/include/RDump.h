/// RDump.h
/// ----------------------------------------------------------------------------
/// rep RDump c++ namespace
/// 23/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------


#ifndef DEF_RDUMP
#define DEF_RDUMP


/// includes
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <string>

namespace RDump {

/// line
/// ----------------------------------------------------------------------------
/// dump a line
void line(int width=80, std::string const& charac="*", 
          std::ostream& flux=std::cout);


} /// end of namespace 


#endif /// RDUMP
