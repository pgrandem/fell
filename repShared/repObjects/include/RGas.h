/// RGas.h
/// ----------------------------------------------------------------------------
/// rep RGas c++ class 
/// Pierre GRANDEMANGE
/// 15/06/2017
/// ----------------------------------------------------------------------------


#ifndef DEF_RGAS
#define DEF_RGAS


/// includes
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
/// root classes
/// rep namespace
#include "RUnits.h"
/// rep classes
#include "RMolecule.h"
#include "RObject.h"

class RFlyer;
class RMachine;
class RBeam;

class RGas : public RObject
{
  /// attribute
  /// **************************************************************************
  private :
  std::vector<RMolecule*> rcompound;      /// [ ]    gas molecular compounds
  std::vector<double>     rfraction;      /// [ ]    fraction of each molecules
  double                  rpressureN2eq;  /// [mbar] N2 equivalent pressure
  double                  rtemperature;   /// [K]    overall gas temperature
  std::string             rgauge;         /// [ ]    gauge used to mesure pS
  
  /// constructor, destructor, copy 
  /// **************************************************************************
  public :
  RGas(std::string const& name);
  RGas(std::string const& name, double pressure, double temperature);
  RGas(RMolecule* const& mol, double pN2eq, double T);
  virtual ~RGas();
  
  
  /// accessors
  /// **************************************************************************
  public :
  std::vector<RMolecule*>   getCompound()     const { return rcompound;     }
  std::string               getGauge()        const { return rgauge;        }
  double                    getPressureN2eq() const { return rpressureN2eq; }
  double                    getTemperature()  const { return rtemperature;  }
  
  /// setGauge
  void setGauge(std:: string const& gaugeModel) { rgauge=gaugeModel; }
  
  /// setPressureN2eq
  /// --------------------------------------------------------------------------
  /// set overall N2 equivalent pressure
  ///   Value generally given by pressure gauge. For each molecule, apply 
  ///   a correction factor given in gauge manual.
  void setPressureN2eq(const double& pN2eq) { rpressureN2eq = pN2eq; };
  
  /// setComposition
  /// --------------------------------------------------------------------------
  ///set molecule composition. Partial pressure indicated in fraction of total.
  void setCompound(std::vector<RMolecule*> const& molecules,
                   std::vector<double>     const& fractions)
    { rcompound = molecules; rfraction = fractions; checkComp();};
                      
  
  
  
  
  /// dump methods
  /// **************************************************************************
  /// dump
  virtual void dump(std::ostream &flux=std::cout) const;
  /// dumpLine
  virtual void dumpLine(std::ostream &flux=std::cout) const;
  
  /// dumpRates (flyer version)
  void dumpRates(RFlyer* const& fly, RMachine* const & mac, 
                 std::ostream &flux=std::cout) const;
  
  /// dumpRates (beam version)
  void dumpRates(RBeam* const& bea, RMachine* const & mac, 
                 std::ostream &flux=std::cout) const;
  
  
  
  
  
  
  
  /// "other" methods
  /// **************************************************************************
  public :
  /// checkComp
  /// --------------------------------------------------------------------------
  /// check if sum of fractional molecular componants = 1;
  bool checkComp() const;
  
  /// density
  /// --------------------------------------------------------------------------
  /// Compute density for each molecule in the gas taking fractionnal
  /// composition and gauge pK correction factor into account.
  std::vector<double> density(bool verbose=true) const;
  
  /// densityFracN2eq
  /// --------------------------------------------------------------------------
  /// Compute N2eq density for each molecule in the gas taking fractionnal
  /// composition into account.
  std::vector<double> densityFracN2eq() const;
  
  /// densityN2eq
  /// --------------------------------------------------------------------------
  /// compute N2eq overall molecule volumic density from p, T
  double densityN2eq() const { return rpressureN2eq/(RUnits::kB*rtemperature); }
  
  /// pPartialN2eq
  /// --------------------------------------------------------------------------
  /// get N2eq partial pressures of molecule constituating the gas
  /// Not physical
  std::vector<double> pPartialN2eq() const;
  
  /// pPartial
  /// --------------------------------------------------------------------------
  /// get partial pressures of molecule constituating the gas
  /// taking into account fractional and pK correction factor for each molecules
  std::vector<double> pPartial(bool verbose=true) const;
  
  
  /// cross sections, interaction rates and gas rates
  /// --------------------------------------------------------------------------
  /// look at RMolecule and RAtom classes...
  /// cross sections per molecules
  std::vector<double> ion_pbarcs(RFlyer* const& fly) const;
  std::vector<double> rs_ebucs (RFlyer* const& fly, RMachine* const& mac) const;
  std::vector<double> rs_lalcs (RFlyer* const& fly, RMachine* const& mac) const;
  std::vector<double> rs_totcs() const;
  /// rates per molecules (per flyer)
  std::vector<double> ion_pbarnu(RFlyer* const& fly) const;
  std::vector<double> rs_ebunu(RFlyer* const& fly, RMachine* const& mac) const;
  std::vector<double> rs_ebunuCarli(RFlyer* const& fly, 
                                    RMachine* const& mac) const;
  std::vector<double> rs_lalnu(RFlyer* const& fly, RMachine* const& mac) const;
  std::vector<double> rs_totnu(RFlyer* const& fly) const;
  /// rates for the gas (per flyer)
  double ion_pbarnuGas(RFlyer* const& fly) const;
  double rs_ebunuGas(RFlyer* const& fly, RMachine* const& mac) const;
  double rs_lalnuGas(RFlyer* const& fly, RMachine* const& mac) const;
  double rs_totnuGas(RFlyer* const& fly) const;
  
  /// rates for the gas per beam (times number of flyers, divided by mac length)
  /// 2017 09 20 : to be change to adapt to new flyer/beam classes
  double ion_pbarnuGas(RBeam* const& beam, RMachine* const& mac) const;
  double rs_ebunuGas(RBeam* const& beam, RMachine* const& mac) const;
  double rs_lalnuGas(RBeam* const& beam, RMachine* const& mac) const;
  double rs_totnuGas(RBeam* const& beam, RMachine* const& mac) const;
  
  /// rsRates()
  /// --------------------------------------------------------------------------
  /// computes beam/gas mean interaction rates and standard deviation(sigma)
  /// sum flyer/gas interaction rates over the beam flyer collection
  /// returns nu [#]/[L]/[T]
  /// two computation loop to get these values
  /// returns a ptr:
  ///   - ptr[0] | ion_pbarnuGas(beam) = sum( ion_pbarnuGas(flyer)  )/macLength
  ///   - ptr[1] | rs_ebunuGas(beam)   = sum( rs_ebunuGas(flyer)    )/macLength
  ///   - ptr[2] | rs_lalnuGas(beam)   = sum( rs_lalnuGas(flyer)    )/macLength
  ///   - ptr[3] | rs_totnuGas(beam)   = sum( rs_totnuGas(flyer)    )/macLength
  ///   - ptr[4] | sigma( ion_pbarnuGas(beam) )
  ///   - ptr[5] | sigma( rs_ebunuGas(beam)   )
  ///   - ptr[6] | sigma( rs_lalnuGas(beam)   )
  ///   - ptr[7] | sigma( rs_totnuGas(beam)   )
  double* rsRates(RBeam* const& beam, RMachine* const& mac) const; 
  
  
  
  
  
};

#endif /// RGAS
