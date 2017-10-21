/// main
/// ----------------------------------------------------------------------------
/// "main" for BIDays program.
/// Pierre GRANDEMANGE
/// 12/06/2017
/// ----------------------------------------------------------------------------

/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <ctime>      /// clock
#include <iostream>
/// root classes
/// rep classes
#include "RDump.h"
/// rep namespaces
/// rep data functions
/// local functions
#include "testFunc.h"
/// namespace
using namespace std;



/// prototypes
/// ****************************************************************************
void hiRep();  /// hello rep test function



/// main program
/// ****************************************************************************
int main()
{
  /// main intro
  /// --------------------------------------------------------------------------
  /// clock main
  const clock_t tstart = clock();   /// main start time
  float tmain(0);                   /// main time counter
  /// intro header
  cout << endl; cout << "fell main() start" << endl; RDump::line(17);
  
  
  /// main : here it goes!
  /// --------------------------------------------------------------------------
  testNewGFB();
  //testBeamDist();
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
  cout << endl; cout << "fell main() end" << endl; RDump::line(35);
  /// main end computation time
  tmain = 1000.*float(clock()-tstart)/CLOCKS_PER_SEC;
  cout << "<main> computation time = " << tmain << "ms" << endl << endl; 
  return 0;
}




/// functions
/// ****************************************************************************

/// hello rep test function
/// ----------------------------------------------------------------------------
/// - dump "hello rep"
void hiRep()
{
  cout << "hello rep" << endl;  
}
