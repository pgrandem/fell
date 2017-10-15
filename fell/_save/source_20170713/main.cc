/// main
/// ----------------------------------------------------------------------------
/// "main" for BIDays program.
/// Pierre GRANDEMANGE
/// 12/06/2017
/// ----------------------------------------------------------------------------

/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <stdlib.h>     /// getenv
/// root classes
/// rep classes
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
  /// intro
  cout << endl;
  cout << "prBIDays main() start" << endl;
  cout << "*********************" << endl;
  
  /// here it goes!
  /// --------------------------------------------------------------------------
  //ecplots();
  //ecoolplot();
  //elenaCycle();
  elenaRate();
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
  
  
  /// end of program
  cout << endl;
  cout << "prBIDays main() end" << endl;
  cout << "*******************" << endl;
  cout << endl;
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
