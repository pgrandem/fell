/// RParticleList.h
/// ----------------------------------------------------------------------------
/// rep RParticleList c++ namespace
/// 13/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------

/// include existing particles ones for all

#ifndef DEF_RPARTICLELIST
#define DEF_RPARTICLELIST


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// rep namespaces
#include "RUnits.h"
/// rep classes
#include "RAtom.h"
#include "RParticle.h"
#include "RMolecule.h"


namespace RParticleList {
/// atoms definition
/// ****************************************************************************
/// list of isotopes by atomic number, mass number
///   - http://physics.nist.gov/cgi-bin/Compositions/
///     stand_alone.pl?ele=&ascii=html&isotype=some
const RAtom* aH1(); 
const RAtom* aC12();
const RAtom* aN14();
const RAtom* aO16();
const RAtom* aNe20();
const RAtom* aAr40();
const RAtom* aKr84();
const RAtom* aXe132();





/// hadrons
/// ****************************************************************************

/// antiproton
/// ----------------------------------------------------------------------------
const RParticle* hpbar();
/// rephadron
/// ----------------------------------------------------------------------------
const RParticle* hrep();






/// leptons
/// ****************************************************************************
/// electron





/// molecules definition (alphabetic order)
/// ****************************************************************************
/// Ar
/// ----------------------------------------------------------------------------
const RMolecule* mAr();

/// CH4
/// ----------------------------------------------------------------------------
const RMolecule* mCH4();
/// CO
/// ----------------------------------------------------------------------------
const RMolecule* mCO();
/// CO2
/// ----------------------------------------------------------------------------
const RMolecule* mCO2();
/// H2
/// ----------------------------------------------------------------------------
const RMolecule* mH2();
/// Kr
/// ----------------------------------------------------------------------------
const RMolecule* mKr();
/// N2
/// ----------------------------------------------------------------------------
const RMolecule* mN2();
/// Ne
/// ----------------------------------------------------------------------------
const RMolecule* mNe();
/// Xe
/// ----------------------------------------------------------------------------
const RMolecule* mXe();


}

#endif /// namespace particleList




