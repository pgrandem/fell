/// localFunctions.cc
/// ----------------------------------------------------------------------------
/// local functions specific to this program folder
/// Pierre Grandemange
/// 09.06.2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
/// root classes
/// rep classes
#include "RAtom.h"
#include "RObject.h"
#include "RParticle.h"
/// rep namespaces
#include "RUnits.h"
#include "RAtoms.h"
/// local functions
#include "testFunc.h"
using namespace std;
using namespace RUnits;


/// category 01
/// ****************************************************************************

/// tAtoms - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RAtoms namespace
void tAtoms()
{
  /// function header
  cout << endl << "tAtoms()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
  cout << "test Atoms namespace with Atom object declaration" << endl;
  RAtom* ato(0);
  
  ato = &RAtoms::hyd; cout << endl << *ato << endl;
  ato = &RAtoms::car; cout << endl << *ato << endl;
  ato = &RAtoms::nit; cout << endl << *ato << endl;
  ato = &RAtoms::oxy; cout << endl << *ato << endl;
  
}


/// testUnit - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RUnit class
void tUnit()
{
  /// No... I don't really know what to do with this class...
  /// UNFINISHED class
}

/// testUnits - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RUnits namespace
void tUnits()
{
  cout << 1.*RUnits::kilogram << endl;
}

/// testAtom - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RAtom class
void tAtom()
{
  /// function header
  cout << "test Atom class function start" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
  
  /// test class, basics
  RAtom* ato1 = new RAtom("repAtom", "Re", 100., 200.);
  cout << endl; ato1->dump();
  
  cout << endl << "test atom get accessors" << endl;
  cout << ato1->getName() << endl;
  cout << ato1->getSymbol() << endl;
  cout << ato1->getMass() << endl;
  cout << ato1->getCharge() << endl;
  cout << ato1->getZ() << endl;
  cout << ato1->getA() << endl;
  
  /// test polymorphism
  cout << endl << "test atom polymorphism (<< operator)" << endl;
  testPol(ato1);
  
  /// test total/geometric cross section
  cout << endl << "test atom geometric cross section" << endl;
  delete ato1; ato1=0;
  ato1 = new RAtom("Nitrogen", "N", 7, 14, 14.007*amu);
  cout << *ato1 << endl;
  
  
  /// test atoms creation in RAtoms.h
  cout << endl << "test atoms creation in RAtoms.h" << endl;
  delete ato1; ato1=0;
  ato1 = &RAtoms::nit;
  cout << *ato1 << endl;
  
 
  
  /// end of function
  cout << endl; cout << "test Atom class function end" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
}


/// testPol - 12.06.2017
/// ----------------------------------------------------------------------------
/// test polymorphisme
void testPol(RObject* const obj)
{
  cout << *obj << endl;
}

/// testRP - 12.06.2017
/// ----------------------------------------------------------------------------
/// test class RParticle
void testRP01()
{
  cout << endl;
  cout << "RObject::nb     " << RObject::nb << endl;
  cout << "RParticle::nb() " << RParticle::nb() << endl;
  cout << endl;
  RParticle::list();
  
  RObject* obj1 = new RObject("obj1");
  RParticle* par1 = new RParticle("p", "proton", 1., 938.e6);
  RParticle* par2 = new RParticle("n", "neutron", 0., 940.e6);
  RParticle* par3 = new RParticle("He2+", "Helium 2+", 2., 4.*940.e6);
  
  cout << endl;
  testPol(obj1);
  cout << endl;
  testPol(par1);
  cout << endl;
  testPol(par2);
  cout << endl;
  testPol(par3);
  
  cout << endl;
  cout << "RObject::nb     " << RObject::nb << endl;
  cout << "RParticle::nb() " << RParticle::nb() << endl;
  cout << endl;
  RParticle::list();
  
  delete par3; par3=0;
  
  cout << endl;
  cout << "RObject::nb     " << RObject::nb << endl;
  cout << "RParticle::nb() " << RParticle::nb() << endl;
  cout << endl;
  RParticle::list();
  
}


/// testRObject - 09.06.2017
/// ----------------------------------------------------------------------------
/// test of object RObject
void testRObject()
{
  /// declare pointers to new objects
  RObject* obj1 = new RObject();
  RObject* obj2 = new RObject();
  
  /// do stuff
  obj1->setName("name1");
  obj2->setName("name2");
  RObject obj3 = *obj1; /// call copy constructor
  RObject obj4 = obj3;  /// call copy constructor
  cout << obj1->getName() << "  " << obj2->getName() << "  " 
       << obj3.getName()  << "  " << obj4.getName()  << "  " 
       << endl;
  
  obj4 = *obj2; /// call "=" operator
  cout << obj1->getName() << "  " << obj2->getName() << "  " 
       << obj3.getName()  << "  " << obj4.getName()  << "  " 
       << endl;
  
  /// release memory
  delete obj1; obj1=0;
  delete obj2; obj2=0;
  
  bool infTest = obj3<obj4; /// test "<" operator
  cout << infTest << endl;
  cout << (obj3<obj4) << endl;
}




