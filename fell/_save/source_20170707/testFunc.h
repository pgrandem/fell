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

/// ecoolplot
/// ----------------------------------------------------------------------------
/// plot ecool graph
void ecoolplot();

/// ecEval
/// ----------------------------------------------------------------------------
/// eval y(x) on electron cooling simulation curve
double ecoolEmit(TGraph* const& tg, double const& x, 
                 double const& yi, double const& yf);

/// ecool(string const& coolingStep)
/// ----------------------------------------------------------------------------
/// make ecooling graph from Gerard Tranquille simulation file (betaCool)
TGraph* ecool(std::string const& coolingStep);


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


/// tRate() - 22/06/2017
/// ----------------------------------------------------------------------------
/// test new atom/mol/gas "old style" for rutherford scattering results
void tRS();

/// tInt() - 20.06.2017
/// ----------------------------------------------------------------------------
/// test Interaction_FLyerGas class
void tInt();

/// tMac() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Machine/Flyer
void tMac();

/// tCSgeom() - 19.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class
/// bilan : on test les cs geom et interaction rate pour comparer avec 
/// resultats 19 Juin 2015 (papier)
//void tCSgeom();  UNFINISHED
// useless without beam class and cs class

/// tgas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class
void tgas();

/// tgas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class
void tgas();

/// molIonCSeval() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Ionisation cross section eval/interpolation function
void molIonCSeval();

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

#endif /// TESTFUNC
