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
#include <ctime>      /// clock_t, clock(), CLOCKS_PER_SEC
#include <iostream>
#include <string>

namespace RDump {

// line
/// ----------------------------------------------------------------------------
/// dump a line
void line(int width=80, std::string const& charac="*", 
          std::ostream& flux=std::cout);


/// timer
/// ----------------------------------------------------------------------------
/// dump time in ms spent since <timerstart>
///   - returns dumped value
float timer(clock_t timerstart, std::ostream &flux=std::cout);

/// timerms
/// ----------------------------------------------------------------------------
/// returns time in ms spent since <clockstart>
float timerms(clock_t timerstart);

/// header
/// ----------------------------------------------------------------------------
/// dump header with a header sized line of <-> characters
/// retruns the title
std::string header(std::string const& header, std::ostream& flux=std::cout);



} /// end of namespace 


#endif /// RDUMP
