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
//#include <stdlib.h>   /// system
#include <iostream>
#include <string>
/// root class
#include "TCanvas.h"

namespace RDump {

/// buildFolder
/// ----------------------------------------------------------------------------
/// build folder <folder> if it does not exist
/// returns true if do not exist, false if already exist
/// https://linux.die.net/man/3/mkdir
bool buildFolder(std::string const& folder);

/// canzenbook3
/// ----------------------------------------------------------------------------
/// return a "root" TCanvas that fills 95% of a quarter of my zenbook3 screen
TCanvas* canzenbook3(std::string const& can="can", 
                     std::string const& title="title");

/// checkFolder
/// ----------------------------------------------------------------------------
/// return 1 if folder exist, 0 otherwise
/// https://www.quora.com/How-do-I-check-if-a-directory-exists-in-Linux
///   -using-C++?share=1
bool checkFolder(std::string const& folder);

/// date
/// ----------------------------------------------------------------------------
/// return a string with date (default=yyyymmdd)
std::string date(std::string const& format="%Y%m%d");

/// header
/// ----------------------------------------------------------------------------
/// dump header with a header sized line of <-> characters
/// retruns the title
std::string header(std::string const& header, std::ostream& flux=std::cout);

/// line
/// ----------------------------------------------------------------------------
/// dump a line
std::string line(int width=80, std::string const& charac="*", 
                 std::ostream& flux=std::cout);

/// pathNow
/// ----------------------------------------------------------------------------
/// returns file path with today's result folder: [...]/results/today/file
std::string pathNow(std::string const& file="", 
                    std::string const& format="%Y%m%d");

/// timer
/// ----------------------------------------------------------------------------
/// dump time in ms spent since <timerstart>
///   - returns dumped value
float timer(clock_t timerstart, std::ostream &flux=std::cout);

/// timerms
/// ----------------------------------------------------------------------------
/// returns time in ms spent since <clockstart>
float timerms(clock_t timerstart);



} /// end of namespace 


#endif /// RDUMP
