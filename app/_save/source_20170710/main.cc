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
//#include <RFlyer.h>
/// rep namespaces
//#include <RUnits.h>
/// rep data functions
#include "elenaCool.h"
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
  cout << endl;
  
  
  /// here it goes!
  /// --------------------------------------------------------------------------
  //ecplots();
  //ecoolplot();
  //elenaCycle();
  //elenaRate();
  //tRS();
  //tInt();
  //tMac();
  //tCSgeom();  /// UNFINISHED
  //tgas();
  //molIonCSeval();
  //molIonisationCS();
  //molList();
  listParticles();
  //tAtom();
  //testRP01();
  //testRObject();
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