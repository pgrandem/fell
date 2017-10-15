/// RParticleList.cc
/// ----------------------------------------------------------------------------
/// rep RParticleList c++ namespace
/// 13/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------

/// include existing particles ones for all

/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
//#include <string>
/// rep classes
#include "RAtom.h"
#include "RMolecule.h"
/// rep namespaces
#include "RParticleList.h"
#include "RUnits.h"


namespace RParticleList {

/// atoms definition
/// ****************************************************************************
/// list of isotopes by atomic number, mass number
///   - http://physics.nist.gov/cgi-bin/Compositions/
///     stand_alone.pl?ele=&ascii=html&isotype=some

/// atomic hydrogen 1
const RAtom* aH1()
{
  const RAtom* atom = new 
        RAtom("Hydrogen", "H", 1, 1, 1.00782503223*RUnits::amu);
  return atom;
}

/// atomic carbon 12
const RAtom* aC12()
{
  const RAtom* atom = new 
        RAtom("Carbon", "C", 6, 12, 12.0000000*RUnits::amu);
  return atom;
}

/// atomic nitrogen 14
const RAtom* aN14()
{
  const RAtom* atom = new 
        RAtom("Nitrogen", "N", 7, 14, 14.00307400443*RUnits::amu);
  return atom;
}

/// atomic oxygen 16
const RAtom* aO16()
{
  const RAtom* atom = new 
        RAtom("Oxygen","O", 8, 16, 15.99491461957*RUnits::amu);
  return atom;
}

/// atomic neon 20
const RAtom* aNe20()
{
  const RAtom* atom = new 
        RAtom("Neon","Ne", 10, 20, 19.9924401762*RUnits::amu);
  return atom;
}

/// atomic argon 40
const RAtom* aAr40()
{
  const RAtom* atom = new 
        RAtom("Argon","Ar", 18, 40, 39.9623831237*RUnits::amu);
  return atom;
}

/// atomic krypton 84
const RAtom* aKr84()
{
  const RAtom* atom = new 
        RAtom("Krypton","Kr", 36, 84, 83.9114977282*RUnits::amu);
  return atom;
}

/// atomic xenon 132
const RAtom* aXe132()
{
  const RAtom* atom = new 
        RAtom("Xenon","Xe", 54, 132, 131.9041550856*RUnits::amu);
  return atom;
}






/// hadrons
/// ****************************************************************************


/// rephadron
/// ----------------------------------------------------------------------------
const RParticle* hrep()
{
  const RParticle* part = new 
        RParticle("rephadron", "repon", -3.1415*RUnits::e, -1000.*RUnits::amu);
  return part;
}

/// antiproton
/// ----------------------------------------------------------------------------
const RParticle* hpbar()
{
  const RParticle* part = new 
        RParticle("antiproton", "pbar", -1.*RUnits::e, 1.0072765*RUnits::amu);
  return part;
}







/// leptons
/// ****************************************************************************

/// electron







/// molecules definition (alphabetic order)
/// ****************************************************************************
/// mass values: 
///   http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440633&Units=SI&Mask=21

/// Ar
/// ----------------------------------------------------------------------------
const RMolecule* mAr()
{
  RMolecule* mol = new RMolecule("argon","Ar","Ar40");
  mol->addAtom((RAtom*)aAr40());
  mol->setMass(39.948*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// CH4
/// ----------------------------------------------------------------------------
const RMolecule* mCH4()
{
  RMolecule* mol = new RMolecule("methane","CH4","C12-H1_4");
  mol->addAtom((RAtom*)aC12()); 
  mol->addAtom((RAtom*)aH1()); 
  mol->addAtom((RAtom*)aH1()); 
  mol->addAtom((RAtom*)aH1()); 
  mol->addAtom((RAtom*)aH1());
  mol->setMass(16.0425*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// CO
/// ----------------------------------------------------------------------------
const RMolecule* mCO()
{
  RMolecule* mol = new RMolecule("carbon_monoxide","CO","C12-O16");
  mol->addAtom((RAtom*)aC12()); 
  mol->addAtom((RAtom*)aO16()); 
  mol->setMass(28.0101*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// CO2
/// ----------------------------------------------------------------------------
const RMolecule* mCO2()
{
  RMolecule* mol = new RMolecule("carbon_dioxide","CO2","C12-O16_2");
  mol->addAtom((RAtom*)aC12()); 
  mol->addAtom((RAtom*)aO16()); 
  mol->addAtom((RAtom*)aO16());
  mol->setMass(44.0095*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// H2
/// ----------------------------------------------------------------------------
const RMolecule* mH2()
{
  RMolecule* mol = new RMolecule("dihydrogen","H2","H1_2");
  mol->addAtom((RAtom*)aH1()); 
  mol->addAtom((RAtom*)aH1()); 
  mol->setMass(2.01588*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// Kr
/// ----------------------------------------------------------------------------
const RMolecule* mKr()
{
  RMolecule* mol = new RMolecule("krypton","Kr","Kr84");
  mol->addAtom((RAtom*)aKr84()); 
  mol->setMass(83.798*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// N2
/// ----------------------------------------------------------------------------
const RMolecule* mN2()
{
  RMolecule* mol = new RMolecule("dinitrogen","N2","N14_2");
  mol->addAtom((RAtom*)aN14()); 
  mol->addAtom((RAtom*)aN14()); 
  mol->setMass(28.0134*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// Ne
/// ----------------------------------------------------------------------------
const RMolecule* mNe()
{
  RMolecule* mol = new RMolecule("neon","Ne","Ne20");
  mol->addAtom((RAtom*)aNe20()); 
  mol->setMass(20.1797*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}

/// Xe
/// ----------------------------------------------------------------------------
const RMolecule* mXe()
{
  RMolecule* mol = new RMolecule("xenon","Xe","Xe132");
  mol->addAtom((RAtom*)aXe132()); 
  mol->setMass(131.293*RUnits::amu);
  mol->setCharge(0.);
  return mol;
}


} /// end of namespace RPArticleList






