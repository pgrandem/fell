/// RMachine.h
/// ----------------------------------------------------------------------------
/// rep RMachine c++ class 
/// Pierre GRANDEMANGE
/// 19/06/2017
/// ----------------------------------------------------------------------------


#ifndef DEF_RMACHINE
#define DEF_RMACHINE


/// includes
/// ----------------------------------------------------------------------------
/// standard library
#include <string>
#include <cmath>
#include <algorithm> /// std::min
/// root classes
/// rep classes
#include "RObject.h"


class RMachine : public RObject
{
  public: 
  /// constructor, destructor, copy 
  /// **************************************************************************
  RMachine();
  RMachine(std::string const& machineName);
  RMachine(std::string const& machineName, double machineLength);
  RMachine(std::string const& name, double length, double accept, double betaT);
  virtual ~RMachine();
  
  
  
  /// accessors
  /// **************************************************************************
  double getAcceptance()  const { return racceptance; 	};
  double getLength()      const { return rlength;     	};
  double getTwissBetaX()  const { return rtwissBeta[0]; };
  double getTwissBetaY()  const { return rtwissBeta[1]; };
  
  
  
  /// methods
  /// **************************************************************************
  double acceptanceXangle() const { return sqrt(racceptance/rtwissBeta[0]); };
  double acceptanceYangle() const { return sqrt(racceptance/rtwissBeta[1]); };
  /// min of acceptanceAngle x/y
  double acceptanceAngle()  const 
    { return std::min(acceptanceXangle(), acceptanceYangle()); };
  
	/// dump
	/// --------------------------------------------------------------------------
	/// dump properties
	void dumpLine(std::ostream &flux=std::cout) const;

  
  /// attribute
  /// **************************************************************************
  protected:
  double racceptance;     /// [length] machine nominal acceptance
  double rlength;         /// [length] machine length
  double rtwissBeta[2];   /// [length] twiss beta value: x and y
  
};

# endif /// RMACHINE

