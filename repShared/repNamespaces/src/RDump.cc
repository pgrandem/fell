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
/// system library
#include <sys/stat.h> /// check folder, mkdir
/// rep classes
/// rep namespaces
#include "RDump.h"
#include "RUnits.h"
//using namespace std;



/// forward declaration for const global variable defined elsewhere
/// ****************************************************************************
extern const std::string resultFolder; /// main timer start



/// namespace
/// ****************************************************************************
namespace RDump {
/// buildFolder
/// ----------------------------------------------------------------------------
/// build folder <folder> if it does not exist
/// returns true if do not exist, false if already exist
/// https://linux.die.net/man/3/mkdir
bool buildFolder(std::string const& folder)
{
  bool cf(checkFolder(folder.c_str()));
  if( !cf ) { 
    std::cout << "RDump::buildfolder -> ";
    std::cout << "building new folder: " << std::endl;
    std::cout << "\"" << folder << "\"" << std::endl;
    mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
  return !cf;
}



/// canzenbook3
/// ----------------------------------------------------------------------------
/// return a "root" TCanvas that fills 95% of a quarter of my zenbook3 screen
TCanvas* canzenbook3(std::string const& name, std::string const& title)
{
  const double cw(864); /// can width  = 0.95 * 1920/2
  const double ch(486); /// can height = 0.95 * 1080/2 
  TCanvas* can = new TCanvas(name.c_str(), title.c_str(), cw, ch);
  can->SetWindowSize( cw+(cw-can->GetWw()), ch+(ch-can->GetWh()) );
  return can;
}



/// checkFolder
/// ----------------------------------------------------------------------------
/// return 1 if folder exist, 0 otherwise
/// https://www.quora.com/How-do-I-check-if-a-directory-exists-in-Linux
///   -using-C++?share=1
bool checkFolder(std::string const& folder)
{
  struct stat statbuf;
  int isDir = 0;
  if( stat(folder.c_str(), &statbuf) != -1 ) {
    if( S_ISDIR(statbuf.st_mode) ) {
      isDir = 1;
    }
  }
  else {
    /* 
    here you might check errno for the reason, ENOENT will mean the path
    was bad, ENOTDIR means part of the path is not a directory, EACCESS     
    would mean you can't access the path. Regardless, from the point of 
    view of your app, the path is not a directory if stat fails. 
    */
  }
  return isDir;
}



/// date
/// ----------------------------------------------------------------------------
/// return a string with date (default=yyyymmdd)
/// http://www.cplusplus.com/reference/ctime/strftime/
std::string date(std::string const& format)
{
  /// variables
  time_t rawtime;       /// time variable
  struct tm * timeinfo; /// calendar...
  char cha[80];         /// output
  /// getting time
  time(&rawtime);       /// get time
  timeinfo = localtime(&rawtime);
  /// build string
  strftime(cha, 80, format.c_str(), timeinfo);
  std::string str = cha;
  /// return
  return str;
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



/// line
/// --------------------------------------------------------------------------
/// dump a line
std::string line(int width, std::string const& charac, std::ostream& flux)
{ 
  std::string str("");
  //for( int i=0; i<width; ++i ) { 
  //  flux << charac; 
  //} flux << '\n'; 
  for( int i=0; i<width; ++i ) { str += charac; } 
  flux << str << '\n';
  return str;
}



/// pathNow
/// ----------------------------------------------------------------------------
/// returns file path with today's result folder: [...]/results/today/file
std::string pathNow(std::string const& file, std::string const& format)
{ 
  const std::string dat(date(format));          /// today
  const std::string now(resultFolder + "/" + dat);
  const bool        bfc(buildFolder(now));
  std::string path(now);
  if( file.size() != 0 ) { path = now + "/" + file; };
  /// debug
  //std::cout << "RDump::pathNow " << "dat " << dat << std::endl;
  //std::cout << "RDump::pathNow " << "now " << now << std::endl;
  //std::cout << "RDump::pathNow " << "bfc " << bfc << std::endl;
  //std::cout << "RDump::pathNow " << "path " << path << std::endl;
  return path;
}



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






} /// end of namespace
