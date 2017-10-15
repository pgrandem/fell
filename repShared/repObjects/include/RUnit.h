/// RUnit.h
/// ----------------------------------------------------------------------------
/// rep RUnit c++ class 
/// Pierre GRANDEMANGE
/// 24/09/2017
/// based on the geant4 model g4UnitSystem
/// ----------------------------------------------------------------------------


#ifndef DEF_RUNIT
#define DEF_RUNIT


/// includes
/// ----------------------------------------------------------------------------
/// standard library
#include <string>
/// root classes
/// rep classes
#include "RObject.h"


class RUnit : public RObject
{
  /// attribute
  /// **************************************************************************
  protected:
  std::string rsymb; /// [ ] unit symbol
  std::string rcate; /// [ ] unit category
  double      rvalu; /// [SI] value in SI unit
 
  /// constructor, destructor, copy 
  /// **************************************************************************
  public: 
  RUnit(std::string name, std::string symbol, std::string category, double val);
  virtual ~RUnit();
  
  /// accessors
  /// **************************************************************************
  std::string getSymbol()   const { return rsymb; }
  std::string getCategory() const { return rcate; }
  double      getValue()    const { return rvalu; }
  
  /// accessors alias
  std::string s() const { return this->getSymbol(); }
  double      v() const { return this->getValue();  }
  /// static methods
  /// **************************************************************************
  //static RUnit tut();
};

double operator*(double val, RUnit const& unit);
double operator*(RUnit const& unit1, RUnit const& unit2);
double operator/(double val, RUnit const& unit);
double operator/(RUnit const& unit1, RUnit const& unit2);

# endif /// RUNIT


