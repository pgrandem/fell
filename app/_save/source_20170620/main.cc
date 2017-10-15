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
/// root classes
/// rep classes
/// rep shared functions
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
  tMac();
  //tCSgeom();  /// UNFINISHED
  //tgas();
  //molIonCSeval();
  //molIonisationCS();
  //molList();
  //tPList();
  //tUnits();
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
