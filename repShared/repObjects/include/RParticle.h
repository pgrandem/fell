/// RParticle.h
/// ----------------------------------------------------------------------------
/// rep RParticle c++ class
/// Pierre Grandemange
/// 12/06/2017
/// ----------------------------------------------------------------------------

#ifndef DEF_RPARTICLE
#define DEF_RPARTICLE


/// includes
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
/// root classes
/// rep classes
#include "RObject.h"



class RParticle : public RObject
{
  public : 
  /// constructor, destructor, copy 
  /// **************************************************************************
  RParticle();
  RParticle(std::string name, std::string symbol);
  RParticle(std::string name, std::string symbol, double charge, double mass);
  virtual ~RParticle();
  
  
  /// accessors
  /// **************************************************************************
  double      getCharge() const   { return rcharge;   };
  double      getMass()   const   { return rmass;     };
  std::string getSymbol() const   { return rsymbol;  };
  
  void  setCharge(double const& charge) { rcharge = charge; };
  void  setMass(double const& mass)     { rmass = mass;     };
  void  setSymbol(double const& symbol) { rsymbol = symbol; };
  
  /// methods
  /// **************************************************************************
  /// dump properties
  /// --------------------------------------------------------------------------
  virtual void dump(std::ostream &flux=std::cout) const;
  void dumpLine(std::ostream &flux=std::cout) const;
  
  /// static methods  
  /// **************************************************************************
  static int  nb(); /// return the number of instance created
  static void list(std::ostream &flux=std::cout); /// return a list of instances
  
  
  /// attributes
  /// **************************************************************************
  protected :
  double        rcharge;  /// [C]   particle charge, use RUnit system 
  double        rmass;    /// [kg]  particle mass
  std::string   rsymbol;  /// [ ]   particle symbol
  
  /// static attributes  
  /// **************************************************************************
  static int                      rnb;   /// number of running instances
  static std::vector<std::string> rlist; /// list of instances (even deleted)
};



/// external operators
/// ****************************************************************************
//std::ostream& operator<<(std::ostream &flux, RParticle const& particle);

#endif /// DEF_RPARTICLE



