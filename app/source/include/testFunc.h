/// localFunctions.h
/// ----------------------------------------------------------------------------
/// local functions specific to this program folder
/// Pierre Grandemange
/// 09.06.2017
/// ----------------------------------------------------------------------------

#ifndef DEF_TESTFUNC
#define DEF_TESTFUNC

/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
/// root classes
#include "TGraph.h"
#include "TH1D.h"
/// rep classes
#include "RObject.h"
/// rep shared functions





/// couple beams, cycles and interactions
/// ****************************************************************************

/// testBeamEmit() - 20171024
/// ----------------------------------------------------------------------------
/// test beam definition with emittance
void testBeamEmit();

/// testBeamDist() - 20171021
/// ----------------------------------------------------------------------------
/// test beam distributions
void testBeamDist();

/// testNewGFB() - 20170921
/// ----------------------------------------------------------------------------
/// test new gas/flyer/beam class 
void testNewGFB();

/// testMeanStd() - 20170921
/// ----------------------------------------------------------------------------
/// test RMath namespace new functions <mean> and <standardDev>
void testRMath();

/// testNewFlyers() - 20170916
/// ----------------------------------------------------------------------------
/// test the new flyer definition, with [x,y,z,px,py,pz]
void testNewFlyers();

/// testBeamProfiles()
/// ----------------------------------------------------------------------------
/// test beam profiles plots
void testBeamProfiles();















/// elena cycle functions
/// ****************************************************************************

/// ecplots()
/// ----------------------------------------------------------------------------
/// elena cycles plots
void ecplots();

/// ecxeg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle x emittance graph (x=horizontal)
TGraph* ecxeg(int const np);

/// ecyeg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle y emittance graph (x=horizontal)
TGraph* ecyeg(int const np);

/// ecekg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle kinetic energy graph 
TGraph* eckeg(int const np);

/// ecmvg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle momentum graph
TGraph* ecmvg(int const np);

/// cg(int ng, int ns, double xs[], double ys[], std::string ymodel);
/// ----------------------------------------------------------------------------
/// cycle graph from array of steps and model to interpol y(x) in between steps
/// ng: nb of graph points
/// ns: nb of steps points
/// xs: x steps array
/// ys: y steps array
/// ymodel: model to interpolate y(x) when xi < x < xi+1
TGraph* cg(int ng, int ns, double xs[], double ys[], std::string const& ymodel);









///beam gas interactions functions
/// ****************************************************************************

/// elenaRate() - 22/06/2017
/// ----------------------------------------------------------------------------
/// compute elena rate with H2
void elenaRate();

/// testRutherfordScattering() - 22.06.2017
/// ----------------------------------------------------------------------------
/// test new atom/mol/gas "old style" for rutherford scattering results
void testRutherfordScattering();

/// testMacFly() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Machine/Flyer classs
void testMacFly();

/// testGas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class
void testGas();

/// testval_csion_pbarGas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Ionisation cross section eval/interpolation function
void testval_csion_pbarGas();

/// plotIon_pbarcs() - 14.06.2017
/// ----------------------------------------------------------------------------
/// plot ion_pbarcs data graphs
/// the folder named "today" must have been created in the results folder 
void plotIon_pbarcs(std::string const& today="toCleanUp");

/// testMolecule() - 14.06.2017
/// ----------------------------------------------------------------------------
/// test RMolecule class
void testMolecule();

/// listParticles() - 07.07.2017
/// ----------------------------------------------------------------------------
/// list particles in particleList namespace
void listParticles();

/// testAtom - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RAtom class
void testAtom();

/// testPol - 12.06.2017
/// ----------------------------------------------------------------------------
/// test polymorphisme with "dump" virtual method
void testPol(RObject* const obj);

/// testParticle() - 12.06.2017
/// ----------------------------------------------------------------------------
/// test class RParticle
void testParticle();

/// testObject - 09.06.2017
/// ----------------------------------------------------------------------------
/// test class RObject
void testObject();

#endif /// TESTFUNC
