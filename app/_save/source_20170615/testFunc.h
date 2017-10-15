/// localFunctions.h
/// ----------------------------------------------------------------------------
/// local functions specific to this program folder
/// Pierre Grandemange
/// 09.06.2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
/// root classes
/// rep classes
#include "RObject.h"
/// rep shared functions


/// category 01
/// ****************************************************************************

/// molIonisationCS() - 14.06.2017
/// ----------------------------------------------------------------------------
/// test Ionisation cross section data retrieving from root file
void molIonisationCS();

/// molList() - 14.06.2017
/// ----------------------------------------------------------------------------
/// list all defined molecules
void molList();

/// tMol() - 14.06.2017
/// ----------------------------------------------------------------------------
/// test RMolecule class
void tMol();

/// tPList - 13.06.2017
/// ----------------------------------------------------------------------------
/// test ParticleList namespace
void tPList();

/// tUnit - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RUnit class
///   No! I don't really know what to do with this class...
///   UNFINISHED class
void tUnit();

/// tUnits - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RUnits namespace
void tUnits();

/// tAtom - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RAtom class
void tAtom();

/// testPol01 - 12.06.2017
/// ----------------------------------------------------------------------------
/// test polymorphisme
void testPol(RObject* const obj);

/// testRP01 - 12.06.2017
/// ----------------------------------------------------------------------------
/// test class RParticle
void testRP01();

/// testObject - 09.06.2017
/// ----------------------------------------------------------------------------
/// test class RObject
void testRObject();


