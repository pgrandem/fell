/// RMolecule.h
/// ----------------------------------------------------------------------------
/// rep RMolecule c++ class 
/// Pierre GRANDEMANGE
/// 12/06/2017
/// ----------------------------------------------------------------------------


#ifndef DEF_RMOLECULE
#define DEF_RMOLECULE


/// includes
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
/// root classes
#include "TGraphErrors.h"
/// rep classes
#include "RParticle.h"
#include "RAtom.h"

class RFlyer;     /// needed for : cross section computation methods
class RMachine;   /// needed for : cross section computation methods

class RMolecule : public RParticle
{
  public:
  
  /// constructors, destructor, copy 
  /// **************************************************************************
  /*
  RMolecule(std::string name, std::string symbol, RAtom* compound, int nbAtoms);
  */
  RMolecule(std::string const& name, std::string const& symbol);
  RMolecule(std::string const& name, std::string const& symbol, 
            std::string const& symbolWithIsotopes);
  RMolecule(std::string const& name, std::string const& symbol, 
            std::string const& symbolWithIsotopes, 
            std::vector<RAtom*> const& comp);
  virtual ~RMolecule();
  
  
  
  /// accessors
  /// **************************************************************************
  void addAtom(RAtom* const& atom) { rcomposition.push_back(atom); };
  
  std::vector<RAtom*> getComposition() const { return rcomposition; };
  std::string getSymbolIsotope() const { return rsymbolIsotope; }
  
  
  
  /// methods
  /// **************************************************************************
  
  /// ion_pbarcs
  /// ----------------------------------------------------------------------------
  /// evaluate pbar cross section for given flyer(ekin)
  /// using ion_pbarcs_data in this class
  ///   -> original graph units: 1.e-20m2 ; keV
  ///   -> input units: RUnits system (SI) -> default=Joule
  ///   -> outut units: RUnits system (SI) -> default=m2
  ///   -> interpolation is linear, ROOT TGraph "eval" method
  double ion_pbarcs(RFlyer* const& flyer) const;
  double ion_pbarcs(double ekin) const;

  /// ion_pbarcs_data
  /// --------------------------------------------------------------------------
  /// retrieve TGraphErrors with pbar ionisation cross section(ekin) data
  /// data from : "J. Phys. B: At. Mol. Opt. Phys. 44 (2011) 122001"
  /// data stored and retrieved in a root file (also in .txt file)
  /// units:  kinE [kEV]
  ///         cs   [1.e-20 m^2]
  TGraphErrors* ion_pbarcs_data() const;
  
  /// rs_ebucs
  /// --------------------------------------------------------------------------
  /// rutherford scattering emittance blow up cross section
  /// IPAC2014, C.Carli, TUPRI028.
  /// sum of atomic cross sections
  double rs_ebucs(RFlyer* const& flyer, RMachine* const& machine) const;
  
  /// rs_lalcs
  /// --------------------------------------------------------------------------
  /// rutherford scattering large angle loss cross section
  /// IPAC2014, C.Carli, TUPRI028.
  double rs_lalcs(RFlyer* const& flyer, RMachine* const& machine) const;
  
  /// rs_totcs
  /// --------------------------------------------------------------------------
  /// rutherford scattering total cross section approximation
  /// IPAC2014, C.Carli, TUPRI028
  /// sum of atomic cross sections
  double rs_totcs(std::string const& model="carli2014") const;
  
  
  
  /// pKfactor
  /// --------------------------------------------------------------------------
  /// get pressure correction factor from pN2eq corresponding to gauge model 
  double pKfactor(double pN2eq, std::string const& gauge="pkr261", 
  								bool verbose=true) const;
  
  /// dump properties
  /// --------------------------------------------------------------------------
  /// dump properties
  virtual void  dump(std::ostream &flux=std::cout) const;
  void          dumpLine(std::ostream &flux=std::cout) const;
  virtual void  dumpcs(RFlyer* fly, RMachine* mac, 
                       std::ostream &flux=std::cout) const;
  
  
  /// attributes
  /// **************************************************************************
  private:
  std::vector<RAtom*> rcomposition; /// [ ] molecule atomic composition
  std::string rsymbolIsotope;             /// [ ] with isotope version of atoms
  
  
  /// static attributes
  /// **************************************************************************
  /// fpbar ionisation cs file name
  static std::string rdataFolder;  
  static std::string rionCSfName;  
};

#endif /// RMOLECULE








