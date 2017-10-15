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
#include "TObject.h"
#include "TGraphErrors.h"
#include "TPad.h"
/// rep classes
#include "RAtom.h"
#include "RMolecule.h"
#include "RObject.h"
#include "RParticle.h"
/// rep namespaces
#include "RUnits.h"
#include "RParticleList.h"
/// local functions
#include "testFunc.h"
using namespace std;
//using namespace RParticleList;
//using namespace RUnits;

/// category 01
/// ****************************************************************************

/// molIonisationCS() - 14.06.2017
/// ----------------------------------------------------------------------------
/// test ionisation cross section data retrieving from root file
void molIonisationCS()
{
  /// function header
  cout << endl; cout << "molIonisationCS()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "check pbar ionisation cs data finding" << endl;
  
  /// do stuff
  RMolecule* mol(0);
  TGraphErrors* tge(0);
  
  mol = RParticleList::mAr();   tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mAr.pdf", "pdf"); 
  gPad->Close();
  
  mol = RParticleList::mCH4();  tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mCH4.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mCO();   tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mCO.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mCO2();  tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mCO2.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mH2();   tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mH2.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mKr();   tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mKr.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mN2();   tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mN2.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mNe();   tge = mol->csIonisationData();
  tge->Draw(); gPad->Print("../results/20170614/cs_mNe.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mXe();   tge = mol->csIonisationData();
  tge->Draw(); 
  gPad->Print("../results/20170614/cs_mXe.pdf", "pdf");
  
  /// end of function
  cout << endl; cout << "molIonisationCS() function end" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}


/// molList() - 14.06.2017
/// ----------------------------------------------------------------------------
/// list all defined molecules
void molList()
{
  /// function header
  cout << endl; cout << "molList()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "dump list of all defined molecules" << endl;
  
  /// dump list of all molecules
  RMolecule* mol(0);
  mol = RParticleList::mAr();   cout << *mol << endl;
  mol = RParticleList::mCH4();  cout << *mol << endl;
  mol = RParticleList::mCO();   cout << *mol << endl;
  mol = RParticleList::mCO2();  cout << *mol << endl;
  mol = RParticleList::mH2();   cout << *mol << endl;
  mol = RParticleList::mKr();   cout << *mol << endl;
  mol = RParticleList::mN2();   cout << *mol << endl;
  mol = RParticleList::mNe();   cout << *mol << endl;
  mol = RParticleList::mXe();   cout << *mol << endl;
  
  /// end of function
  cout << endl; cout << "molList() function end" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
}

/// tMol() - 14.06.2017
/// ----------------------------------------------------------------------------
/// test RMolecule class
void tMol()
{
  /// function header
  cout << endl << "tMol()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
  cout << "test Molecule class \"tMol\" function start" << endl;
  
  /// do stuff
  /// the following is not working
  /*RAtom ato[3];
  ato[0] = RParticleList::aC12;
  ato[1] = RParticleList::aO16;
  ato[2] = RParticleList::aO16;
  RAtom* apt(&ato[0]);
  RMolecule* mol = new RMolecule("Carbon dioxide", "CO2", apt, 3);
  vector<RAtom*> vpt = mol->getComposition();
  for( int i=0; i<vpt.size(); ++i ) { cout << *vpt[i] << endl; }
  */
  
  /// syntax 1
  /// for putting a vector of "const RAtom" insted of "RAtom"
  //RAtom* apt = new RAtom;
  //vector<RAtom*> vap; *apt=*RParticleList::aC12; vap.push_back(apt);
  
  /// syntax 2
  //vap.push_back((RAtom*)RParticleList::aO16);
  //vap.push_back(RParticleList::aO16);
  //for( int i=0; i<vap.size(); ++i ) { cout << *vap[i] << endl; }
  
  
  /// syntax 3
  /*RMolecule* mol = new RMolecule("repMol", "rM", "r30M0");
  mol->addAtom((RAtom*)RParticleList::aN14);
  mol->addAtom((RAtom*)RParticleList::aH1);
  mol->addAtom((RAtom*)RParticleList::aH1);
  mol->addAtom((RAtom*)RParticleList::aH1);
  mol->dump();
  */
  
  
  
  /// end of function
  cout << endl; cout << "test Molecule class \"tMol\" function end" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
}

/// tPList - 13.06.2017
/// ----------------------------------------------------------------------------
/// test ParticleList namespace
/*void tPList()
{
  /// function header
  cout << endl << "tPList()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
  cout << "test ParticleListe namespace with Atom object declaration" << endl;
  
  RAtom* ato(0);
  ato = RParticleList::aH1;    cout << endl << *ato << endl;
  ato = RParticleList::aC12;   cout << endl << *ato << endl;
  ato = RParticleList::aN14;   cout << endl << *ato << endl;
  ato = RParticleList::aO16;   cout << endl << *ato << endl;
  ato = RParticleList::aNe20;  cout << endl << *ato << endl;
  ato = RParticleList::aAr40;  cout << endl << *ato << endl;
  ato = RParticleList::aKr84;  cout << endl << *ato << endl;
  ato = RParticleList::aXe132; cout << endl << *ato << endl;
  
}*/


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
  ato1 = new RAtom("Nitrogen", "N", 7, 14, 14.007*RUnits::amu);
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




