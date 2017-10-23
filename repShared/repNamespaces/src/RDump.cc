/// RDump.cc
/// ----------------------------------------------------------------------------
/// rep RDump c++ namespace
/// 23/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <ctime>      /// clock_t, clock(), CLOCKS_PER_SEC
#include <iostream>
#include <string>
/// rep classes
/// rep namespaces
#include "RDump.h"
#include "RUnits.h"
//using namespace std;


namespace RDump {

/// line
/// --------------------------------------------------------------------------
/// dump a line
void line(int width, std::string const& charac, std::ostream& flux)
{ for( int i=0; i<width; ++i ) { flux << charac; } flux << '\n'; }



/// timer
/// ----------------------------------------------------------------------------
/// dump time in ms spent since <timerstart>
/// calls RDump::timerms()
///   - returns dumped value
float timer(clock_t clockstart, std::ostream &flux)
{
  /// dumped value and output
  auto timer(timerms(clockstart));  
  /// save previous stream format settings
  std::ios::fmtflags cf(flux.flags());      /// save flux current flags
  std::streamsize    cp(flux.precision());  /// save flux current precision
  /// set new stream format settings
  flux.flags(std::ios::left | std::ios::fixed); /// left, fixed notation
  flux.precision(6);
  
  /// do the stuff
  flux << "gtimer=" << timer << RUnits::ms.s() << '\n';
  
  /// restore previous stream format settings
  flux.flags(cf);        /// restore flags
  flux.precision(cp);    /// restore precision
  /// returns stuff
  return timer;
}



/// timerms
/// ----------------------------------------------------------------------------
/// returns time in ms spent since <clockstart>
float timerms(clock_t clockstart)
{
  return 1000.*float(clock()-clockstart)/CLOCKS_PER_SEC;
}



/// header
/// ----------------------------------------------------------------------------
/// dump header with a header sized line of <-> characters
/// retruns <header>
std::string header(std::string const& header, std::ostream& flux)
{
  flux << '\n' << header << '\n';
  line(header.size(), "-", flux);
  return header;
}


} /// end of namespace
