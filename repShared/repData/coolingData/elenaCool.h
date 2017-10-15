/// elenaCool.h
/// ----------------------------------------------------------------------------
/// to deal with cooling simulation files from Gerard Tranquille (betaCool)
/// Pierre Grandemange
/// 07.07.2017
/// ----------------------------------------------------------------------------


#ifndef DEF_ELENACOOL
#define DEF_ELENACOOL

/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
/// root classes
#include "TGraph.h"
/// rep classes
/// rep namespaces



/// ecool(string const& coolingStep)
/// ----------------------------------------------------------------------------
/// get ecooling graph from Gerard Tranquille simulation file (betaCool)
TGraph* ecool(std::string const& coolingStep);

/// scool(string const& coolingStep, double xi, double yi double yf)
/// ----------------------------------------------------------------------------
/// get scaled graph from ecool - scaled ecool - scool
/// scale curve to input param zi and zf, and translate to xi
TGraph* scool(std::string const& coolingStep, double xi, double yi, double yf);

/// vcool
/// ----------------------------------------------------------------------------
/// eval y(x) on electron cooling scaled curve
double vcool(TGraph* const& scool, double const& x);

/// ecoolplot
/// ----------------------------------------------------------------------------
/// plot ecool graph
void ecoolplot(std::string pathOut);



#endif /// ELENACOOL




