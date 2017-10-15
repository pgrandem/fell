/// RAtom.h
/// ----------------------------------------------------------------------------
/// rep RAtom c++ class 
/// Pierre GRANDEMANGE
/// 12/06/2017
/// ----------------------------------------------------------------------------


#ifndef DEF_RATOM
#define DEF_RATOM


/// includes
/// ----------------------------------------------------------------------------
/// standard library
/// root classes
/// rep classes
#include "RParticle.h"


class RFlyer;
class RMachine;

class RAtom : public RParticle
{
  public :
  /// constructor, destructor, copy 
  /// **************************************************************************
  RAtom();
  RAtom(std::string name, std::string symbol);
  RAtom(std::string name, std::string symbol, int Z , int A);
  RAtom(std::string name, std::string symbol, int Z, int A, double mass);
  RAtom(std::string name, std::string symbol, int Z, int A, 
        double mass, double charge);
  virtual ~RAtom();
  
  
  
  /// accessors
  /// **************************************************************************
  double  getA() const   { return rA;   };
  double  getZ() const   { return rZ;   };
  
  void setA(const double& A) { rA = A; };
  void setZ(const double& Z) { rZ = Z; };
  
  
  
  /// methods
  /// **************************************************************************
  
  /// __________________________________________________________________________
  /// rs : rutherford scattering methods.
  /// -> IPAC2014, TUPRI028, Carli modified by me
  
  /// rs_ebucl()
  /// --------------------------------------------------------------------------
  /// rutherford scattering emittance blow up coulomb logarithm
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_ebucl(RFlyer* const& flyer, RMachine* const& machine) const;
  
  /// rs_ebucs()
  /// --------------------------------------------------------------------------
  /// rutherford scattering emittance blow up cross section
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_ebucs(RFlyer* const& flyer, RMachine* const& machine) const;
  
  /// rs_ebuip()
  /// --------------------------------------------------------------------------
  /// rutherford scattering emittance blow up impact parameter
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_ebuip(RFlyer* const& flyer, RMachine* const& machine) const;
  
  /// rs_lalcs()
  /// --------------------------------------------------------------------------
  /// rutherford scattering large angle loss cross section
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_lalcs(RFlyer* const& flyer, RMachine* const& machine) const;
  
  /// rs_lalip()
  /// --------------------------------------------------------------------------
  /// rutherford scattering large angle loss impact parameter
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_lalip(RFlyer* const& flyer, RMachine* const& machine) const;
  
  /// rs_totcs
  /// --------------------------------------------------------------------------
  /// rutherford scattering total cross section approximation
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_totcs(std::string const& model="carli2014") const;
  
  /// rs_totra
  /// --------------------------------------------------------------------------
  /// rutherford scattering total cross section associated radius approximation
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_totra(std::string const& model="carli2014") const;
  
  /// __________________________________________________________________________
  
  
  
  
  /// dump
  /// --------------------------------------------------------------------------
  /// dump properties
  virtual void dump(std::ostream &flux=std::cout) const;
  virtual void dumpLine(std::ostream &flux=std::cout) const;
  virtual void dumpcs(RFlyer* fly, RMachine* mac, 
                      std::ostream &flux=std::cout) const;
  
  
  /// static methods  
  /// **************************************************************************
  
  /// rsip
  /// --------------------------------------------------------------------------
  /// rutherford scattering impact parameter for given angle phi
  /// -> IPAC2014, TUPRI028, Carli
  static double rsip(double Z, double z, double m, double v, double phi);
  
  
  
  /// attributes
  /// **************************************************************************
  private :
  int    rA;       /// [ ]          atomic mass number
  int    rZ;       /// [ ]          atomic number
};

#endif 
