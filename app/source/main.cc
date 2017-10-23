/// main
/// ----------------------------------------------------------------------------
/// "main" for BIDays program.
/// Pierre GRANDEMANGE
/// 12/06/2017
/// ----------------------------------------------------------------------------

/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <ctime>      /// clock_t, clock(), CLOCKS_PER_SEC
#include <iostream>
//#include <stdlib.h>   /// getenv
/// root classes
/// rep classes
#include "RDump.h"
/// rep namespaces
#include "RUnits.h"
/// rep data functions
/// local functions
#include "testFunc.h"
/// namespace
using namespace std;


/// define initialized const global variable (external linkage)
/// ****************************************************************************
/// main timer
extern const clock_t gtimerstart(clock());
/// results folder
extern const std::string resultFolder("/home/rep/programming/fell/app/results");
/// data folder


/// prototypes
/// ****************************************************************************
void hiRep();  /// hello rep test function


/// main program
/// ****************************************************************************
int main()
{
  /// main intro
  /// --------------------------------------------------------------------------
  /// intro header
  cout << endl; cout << "fell main() start" << endl; RDump::line(17);
  
  
  /// main : here it goes!
  /// --------------------------------------------------------------------------
  testBeamDist();
  //testNewGFB();
  //testRMath();
  //testNewFlyers();
  //testBeamProfiles();
  //ecplots();
  //ecoolplot();	/// from data folder
  //elenaRate();
  //testRutherfordScattering();
	//testMacFly();
  //testGas();
  //testval_csion_pbarGas();
  //plotIon_pbarcs(/*today*/"20170712");
	//testMolecule();
  //listParticles();
  //testAtom();
  //testParticle();
  //testObject();
  //hiRep();
  /// --------------------------------------------------------------------------
  
  
  /// main outro
  /// --------------------------------------------------------------------------
  /// outro title
  cout << endl; cout << "fell main() end" << endl; RDump::line(15);
  /// main end computation time
  RDump::timer(gtimerstart);
  /// main return
  return 0;
}




/// functions
/// ****************************************************************************

/// hiRep()
/// ----------------------------------------------------------------------------
/// hello rep test function
/// - dump "hello rep"
void hiRep()
{
  cout << "hello rep" << endl;  
}
