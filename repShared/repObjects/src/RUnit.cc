/// RUnit.cc
/// ----------------------------------------------------------------------------
/// rep RUnit c++ class 
/// Pierre GRANDEMANGE
/// 24/09/2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
/// root classes
/// rep namespaces
//#include "RUnits.h"
/// rep classes
//#include "RObject.h"
#include "RUnit.h"
using namespace std;


/// constructors, destructor, copy 
/// ****************************************************************************
RUnit::RUnit(string name, string symbol, string category, double value) :
  RObject(name), rsymb(symbol), rcate(category), rvalu(value)
{}

RUnit::~RUnit()
{}





/// operators, external
/// ****************************************************************************
double operator*(double val, RUnit const& unit) { return val*unit.getValue(); }
double operator*(RUnit const& unit1, RUnit const& unit2)
{ 
  return unit1.getValue()*unit2.getValue();
}

double operator/(double val, RUnit const& unit) { return val/unit.getValue(); }
double operator/(RUnit const& unit1, RUnit const& unit2)
{
  return unit1.getValue()/unit2.getValue();
}
