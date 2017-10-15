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
//#include "RFlyer.h"
#include "RGas.h"
#include "RInteraction_FlyerGas.h"
#include "RMachine.h"
#include "RMolecule.h"
#include "RObject.h"
#include "RParticle.h"
/// rep namespaces
#include "RDump.h"
#include "RParticleList.h"
#include "RUnits.h"
/// local functions
#include "testFunc.h"
using namespace std;
//using namespace RParticleList;
//using namespace RUnits;

/// category 01
/// ****************************************************************************

/// elenaRate() - 22/06/2017
/// ----------------------------------------------------------------------------
/// compute elena rate with H2
void elenaRate()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl; cout << "elenaRate() - 22.06.2017" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "elena rate with beam" << endl;
  
  
  /// beam parameters
  double npb=1e7;
  
  /// gas
  RGas* gas = new RGas(RParticleList::mH2(), 4.e-12*RUnits::mbar, 300.);
  
  /// flyer
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar, 0.037);
  
  /// machine
  RMachine* mac = new RMachine("elena", 75.*RUnits::um, 3., 30.4055);
  
  /// dump ionisation/cm/s
  double ion_pbar_s = gas->ion_pbarnuGas(fly);
  double ion_s = npb * ion_pbar_s;
  
  cout << ion_s << endl;
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "elenaRate() function end " << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}





/// tRS() - 22.06.2017
/// test new atom/mol/gas "old style" for rutherford scattering results
void tRS()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl; cout << "tRS() - 22.06.2017" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "test new atom/mol/gas for rutherford scattering results" << endl;
  
  
  
  /// do stuff
  /// --------------------------------------------------------------------------
  
  /// define a gas
  double pre(4.e-12*RUnits::mbar);  /// overall pressure
  double tem(300.*RUnits::K);       /// overall temperature
  vector<RMolecule*> mov(0);        /// molecules
  vector<double> mfv;               /// composition fract of molecules
  //mov.push_back(RParticleList::mH2());  mfv.push_back(0.50);
  //mov.push_back(RParticleList::mN2());  mfv.push_back(0.50);
  mov.push_back(RParticleList::mN2());  mfv.push_back(1.00);
  RGas* gas = new RGas("test gas", pre, tem);
  gas->setCompound(mov, mfv);
  /// dump the gas
  //cout << endl; gas->dump(cout, "pkr261");
  
  /// define machine
  RMachine* mac = new 
    RMachine("elena", 75.*RUnits::um, 3.*RUnits::m, 30.4055*RUnits::m);
  //cout << endl;
  //cout << "elena mean acceptance angle: ";
  //cout << mac->acceptanceAngle() << endl;
  
  /// define flyers
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar, 0.037);
  //cout << endl; fly->dump();
  
  
  /// atom dumpcs test
  RAtom* ato = gas->getCompound()[0]->getComposition()[0];
  //cout << endl; ato->dumpcs(fly, mac);
  
  
  /// atom dumpcs test
  //for( int i=0; i<gas->getCompound().size(); ++i ) {
  //  cout << endl;
  //  gas->getCompound()[i]->dumpcs(fly, mac);
  //}
  
  /// dump rates
  cout << endl;
  gas->dumpRates(fly, mac);
  
  delete gas;
  RGas* ggg = new RGas(RParticleList::mH2(), pre, tem);
  gas = ggg;
  
  /// dump rates
  cout << endl;
  gas->dumpRates(fly, mac);
  
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "tRS() function end " << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}








/// tInt() - 20.06.2017
/// ----------------------------------------------------------------------------
/// test Interaction_FLyerGas class
void tInt()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl; cout << "tInt() - 20.06.2017" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "test Interaction_FLyerGas class" << endl;
  
  
  
  /// do stuff
  /// --------------------------------------------------------------------------
  /// define flyers
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar, 0.037);
  cout << endl; fly->dump();
  
  /// define gases
  double pre = 4.e-12*RUnits::mbar;
  double tem = 300.*RUnits::K;
  RGas* dng = new RGas(RParticleList::mN2(), pre, tem); /// dinitrogen gas
  RGas* dhg = new RGas(RParticleList::mH2(), pre, tem); /// dihydrogen gas
  cout << endl; dng->dump(cout, "pkr261"); dhg->dump(cout, "pkr261");
  
  /// the following syntax as been abandonned by including ruth scat cs in
  /// the atom class.
  
  /// define interactions
  /*RInteraction_FlyerGas* npi = new RInteraction_FlyerGas("npi", fly, dng);
  RInteraction_FlyerGas* hpi = new RInteraction_FlyerGas("hpi", fly, dhg);
  //cout << endl; npi->dump();
  //cout << endl; hpi->dump();
  
  /// cross section
  double lossAngle(1.e-3);
  
  cout << endl ;
  cout << npi->getGas()->getName() << " & ";
  cout << npi->getFlyer()->getSymbol() << " @ beta=";
  cout << npi->getFlyer()->getBeta() << " " << endl;
  cout << "sca ="  << npi->sca_csGas() << " ";
  cout << "lasl =" << npi->lasl_csGas(lossAngle) << " ";
  cout << "ebu ="  << npi->ebu_csGas(lossAngle) << " | ";
  cout << "ion ="  << npi->ion_csGas() << " ";
  cout << endl;
  vector<double> laslipN  = npi->lasl_ip(lossAngle);
  vector<double> laslcsN  = npi->lasl_cs(lossAngle);
  vector<double> ebuclN   = npi->ebu_cl(lossAngle);
  vector<double> ebuipN   = npi->ebu_ip(lossAngle);
  vector<double> ebucsN   = npi->ebu_cs(lossAngle);
  vector<double> totcsN   = npi->tot_cs();
  for( int i=0; i<laslipN.size(); ++i ) { 
    cout << "laslipN[i]" << " " << laslipN[i] << "   ";
    cout << "laslcsN[i]" << " " << laslcsN[i] << "   " << endl;
    cout << "ebuclN[i]" << " " << ebuclN[i] << "   ";
    cout << "ebuipN[i]" << " " << ebuipN[i] << "   ";
    cout << "ebucsN[i]" << " " << ebucsN[i] << "   " << endl;
    cout << "totcsN[i]" << " " << totcsN[i] << "   ";
    cout << endl;
  }
  cout << endl ;
  cout << npi->getGas()->getName() << " & ";
  cout << npi->getFlyer()->getSymbol() << " @ beta=";
  cout << npi->getFlyer()->getBeta() << " " << endl;
  cout << "sca_csGas="  << hpi->sca_csGas() << " ";
  cout << "lasl_csGas=" << hpi->lasl_csGas(lossAngle) << " ";
  cout << "ebu_csGas="  << hpi->ebu_csGas(lossAngle) << " | ";
  cout << "ion_csGas="  << hpi->ion_csGas() << " ";
  cout << endl;
  vector<double> laslipH  = hpi->lasl_ip(lossAngle);
  vector<double> laslcsH  = hpi->lasl_cs(lossAngle);
  vector<double> ebuclH   = hpi->ebu_cl(lossAngle);
  vector<double> ebuipH   = hpi->ebu_ip(lossAngle);
  vector<double> ebucsH   = hpi->ebu_cs(lossAngle);
  vector<double> totcsH   = hpi->tot_cs();
  for( int i=0; i<laslipH.size(); ++i ) { 
    cout << "laslipH[i]" << " " << laslipH[i] << "   ";
    cout << "laslcsH[i]" << " " << laslcsH[i] << "   " << endl;
    cout << "ebuclH[i]" << " " << ebuclH[i] << "   ";
    cout << "ebuipH[i]" << " " << ebuipH[i] << "   ";
    cout << "ebucsH[i]" << " " << ebucsH[i] << "   " << endl;
    cout << "totcsH[i]" << " " << totcsH[i] << "   ";
    cout << endl;
  }
  */
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "tInt() function end " << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}

/// tMac() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Machine/Flyer class
void tMac()
{
  /// function header
  cout << endl; cout << "tMac()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "test Machine/Flyer class" << endl;
  
  
  
  /// middle
  RFlyer* fly = new RFlyer("antiproton", "pbar", 
                           -1.*RUnits::e, 1.0072765*RUnits::amu, 0.);
  cout << endl;
  cout << "symbol : " << fly->getSymbol() << endl;
  cout << "name   : " << fly->getName() << endl;
  cout << "mass   : " << fly->getMass()/RUnits::MeV_c2  << " MeV/c2"  << endl;
  cout << "e0     : " << fly->e0()/RUnits::MeV          << " MeV"     << endl;
  cout << "charge : " << fly->getCharge()/RUnits::e     << " e"       << endl;
  
  fly->applyEkin(5.3*RUnits::MeV);
  cout << endl;
  cout << "beta     : " << fly->getBeta()                 << endl;
  cout << "gamma    : " << fly->gamma()                   << endl;
  cout << "ekin     : " << fly->ekin()/RUnits::keV        << " kev"   << endl;
  cout << "momentum : " << fly->momentum()/RUnits::MeV_c  << " MeV/c" << endl;
  
  fly->applyMomentum(35.*RUnits::MeV_c);
  cout << endl;
  cout << "beta     : " << fly->getBeta()                 << endl;
  cout << "gamma    : " << fly->gamma()                   << endl;
  cout << "ekin     : " << fly->ekin()/RUnits::keV        << " kev"   << endl;
  cout << "momentum : " << fly->momentum()/RUnits::MeV_c  << " MeV/c" << endl;
  
  fly->applyEkin(100.*RUnits::keV);
  cout << endl;
  cout << "beta     " << fly->getBeta()                 << endl;
  cout << "gamma    " << fly->gamma()                   << endl;
  cout << "ekin     " << fly->ekin()/RUnits::keV        << " kev"   << endl;
  cout << "momentum " << fly->momentum()/RUnits::MeV_c  << " MeV/c" << endl;
  
  
  
  /// end of function
  cout << endl; cout << "tMac() function end " << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}


/// tCSgeom() - 19.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class
/// bilan : on test les cs geom et interaction rate pour comparer avec 
/// resultats 19 Juin 2015 (papier)
/*void tCSgeom()
{
  /// function header
  cout << endl; cout << "tCSgeom()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "test gas geom cross section to compare old computatioons" << endl;
  
  
  /// pressure and temperature
  const double pre = 4.e-12*RUnits::mbar;
  const double tem = 300.*RUnits::K;
  
  /// atomic nitrogen gas
  vector<RAtom*>  nia(1); nia[0] = (RAtom*)RParticleList::aN14; /// ni atom
  RMolecule* nim = new RMolecule("nitrogen", "N", "N14", nia);  /// ni molec
  RGas* nig = new RGas(nim, pre, tem);  /// ni gas
  nig->dump(cout, "pkr261");
  
  /// dinitrogen gas
  RGas* dng = new RGas(RParticleList::mN2(), pre, tem); /// dinitrogen gas
  dng->dump(cout, "pkr261");
  
  /// dihydrogen gas
  RGas* dhg = new RGas(RParticleList::mN2(), pre, tem); /// dihydrogen gas
  dhg->dump(cout, "pkr261");
  
  
  /// end of function
  cout << endl; cout << "tCSgeom() function end " << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}*/


/// tgas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class
void tgas()
{
  /// function header
  cout << endl; cout << "tgas()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "test Gas class" << endl;
  
  /// gas properties
  double pre(4.e-12*RUnits::mbar);  /// overall pressure
  double tem(298.*RUnits::K);       /// overall temperature
  
  /// define a gas
  vector<RMolecule*> mov(0);        /// molecules
  vector<double> mfv;               /// composition fract of molecules
  mov.push_back(RParticleList::mH2());  mfv.push_back(0.91);
  mov.push_back(RParticleList::mCO2()); mfv.push_back(0.05);
  mov.push_back(RParticleList::mCO());  mfv.push_back(0.03);
  mov.push_back(RParticleList::mCH4()); mfv.push_back(0.01);
  RGas* gas = new RGas("elenaVac", pre, tem);
  gas->setCompound(mov, mfv);
  
  /// dump the gas
  cout << endl << "test elena gas with molecules mixing" << endl;
  gas->dump(cout, "pkr261");
  
  /// simple H2 gas
  //RMolecule* mol = RParticleList::mH2();
  //RGas* sim = new RGas(mol, pre, tem);
  RGas* sim = new RGas(RParticleList::mH2(), pre, tem);
  cout << endl << "test simple H2 gas " << endl;
  sim->dump(cout, "pkr261");
  
  
  
  /// end of function
  cout << endl; cout << "test Gas class" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}


/// molIonCSeval() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Ionisation cross section eval/interpolation function
void molIonCSeval()
{
  /// function header
  cout << endl; cout << "molIonCSeval()" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
  cout << "check pbar ionisation cs evaluation" << endl;
  
  /// do stuff 
  cout << (RParticleList::mCO2())->ion_pbarcs(500.*RUnits::keV) 
       << "m2" << endl;
  
  /// end of function
  cout << endl; cout << "molIonCSeval() function end" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}

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
  
  mol = RParticleList::mAr();   tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mAr.pdf", "pdf"); 
  gPad->Close();
  
  mol = RParticleList::mCH4();  tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mCH4.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mCO();   tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mCO.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mCO2();  tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mCO2.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mH2();   tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mH2.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mKr();   tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mKr.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mN2();   tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mN2.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mNe();   tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mNe.pdf", "pdf");
  gPad->Close();
  
  mol = RParticleList::mXe();   tge = mol->ion_pbarcs_data();
  tge->Draw(); gPad->Print("../results/20170614/cs_mXe.pdf", "pdf");
  gPad->Close();
  
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




