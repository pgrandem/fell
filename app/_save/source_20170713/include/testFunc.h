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
/// rep classes
#include "RObject.h"
/// rep shared functions


/// category 01
/// ****************************************************************************

/// ecplots()
/// ----------------------------------------------------------------------------
/// elena cycle plots
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




/// sequ() - 05/07/2017
/// ----------------------------------------------------------------------------
/// function used in elenaCycle
/// - tsi : time sequence init
/// - tsf : time sequence final
/// - tsv : time step value
/// - vsi : value sequence init
/// - vsf : value sequence final
/// - vsb : value sequence before
/// - bsc : boolean sequence checker
/// - return value evaluated at this step in sequence
double sequ(double tsi, double tsf, double tsv,
               double vsi, double vsf, double vsb, bool& bsc);





/// elenaCycle() - 30/06/2017
/// ----------------------------------------------------------------------------
/// compute elena rate during cycle
void elenaCycle();




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
