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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
/// rep classes
#include "RAtom.h"
#include "RFlyer.h"
#include "RGas.h"
#include "RMachine.h"
#include "RMolecule.h"
#include "RObject.h"
#include "RParticle.h"
/// rep namespaces
#include "RDump.h"
#include "RParticleList.h"
#include "RUnits.h"
#include "RMath.h"
/// local functions
#include "testFunc.h"
using namespace std;
//using namespace RParticleList;
//using namespace RUnits;

/// category 01
/// ****************************************************************************

/// sequ() - 05/07/2017
/// ----------------------------------------------------------------------------
/// sequencer, function used in elenaCycle
/// - tsi : time sequence init
/// - tsf : time sequence final
/// - tsv : time step value
/// - vsi : value sequence init
/// - vsf : value sequence final
/// - vsb : value sequence before
/// - bsc : boolean sequence checker
/// - return value evaluated at this step in sequence
double sequ(double tsi, double tsf, double tsv,
               double vsi, double vsf, double vsb, 
               bool& bsc)
{
  double val(-1);
  double ramp( (vsf-vsi)/(tsf-tsi) );
  double init( vsb );
  if( bsc==false ) { init = vsi; bsc=true; cout << bsc << endl;}
  val = init  + ramp*tsv;
  return val;
}

/// elenaCycle() - 30/06/2017
/// ----------------------------------------------------------------------------
/// ionistations/m/s as function of time
void elenaCycle()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;   cout << "elenaCycle() - 30/06/2017" << endl;
  RDump::line();  cout << "ionistations/m/s as function of time" << endl;
  
  
  /// do the stuff
  /// --------------------------------------------------------------------------
  /// function const param
  const int ns(27);     /// number of samples
  
  /// cycle parameters
  /// ----------------
  
  /// cycle steps : values of interest
  const int cns(10);  /// number of steps in cycle
  double tcs[cns];    /// time cycle steps array
  double kes[cns];    /// kinetic energy steps array
  double mvs[cns];    /// momentum steps array
  double exs[cns];    /// horizontal emittance steps array
  double eys[cns];    /// vertical emittance steps array
  bool   scs[cns];    /// sequence checker steps array
  
  /// time steps - tc = time cycle
  const double tc_inj(0.*RUnits::s);  /// injection;
  const double tc_fdb(1.*RUnits::s);  /// first deceleration begin
  const double tc_fde(5.*RUnits::s);  /// first deceleration end
  const double tc_feb(6.*RUnits::s);  /// first electron cooling start
  const double tc_fee(16.*RUnits::s); /// first electron cooling end
  const double tc_sdb(17.*RUnits::s); /// second deceleration begin
  const double tc_sde(21.*RUnits::s); /// second deceleration end
  const double tc_seb(22.*RUnits::s); /// second electron cooling start
  const double tc_see(25.*RUnits::s); /// second electron cooling end
  const double tc_ext(26.*RUnits::s); /// extraction (rebunched);
  
  /// energy steps - ke = kinetic energy
  const double ke_inj(5.3*RUnits::MeV);   /// injection
  const double ke_fdb(5.3*RUnits::MeV);   /// first deceleration begin
  const double ke_fde(653.*RUnits::keV);  /// first deceleration end
  const double ke_feb(653.*RUnits::keV);  /// first electron cooling start
  const double ke_fee(653.*RUnits::keV);  /// first electron cooling end
  const double ke_sdb(653.*RUnits::keV);  /// second deceleration begin
  const double ke_sde(100.*RUnits::keV);  /// second deceleration end
  const double ke_seb(100.*RUnits::keV);  /// second electron cooling start
  const double ke_see(100.*RUnits::keV);  /// second electron cooling end
  const double ke_ext(100.*RUnits::keV);  /// extraction (rebunched);
  
  /// momentum steps - mv = momentum (mass*velocity)
  /// should correspond to magnetic cycle steps:
  /// -> B*rho # momentum (# = proportional)
  /// shoud be a linear dependance of time
  const double m0 = 1.0072765*RUnits::amu;
  const double mv_inj( RFlyer::ekinToMomentum(ke_inj, m0) );   
  const double mv_fdb( RFlyer::ekinToMomentum(ke_fdb, m0) );  
  const double mv_fde( RFlyer::ekinToMomentum(ke_fde, m0) );
  const double mv_feb( RFlyer::ekinToMomentum(ke_feb, m0) );  
  const double mv_fee( RFlyer::ekinToMomentum(ke_fee, m0) );
  const double mv_sdb( RFlyer::ekinToMomentum(ke_sdb, m0) );
  const double mv_sde( RFlyer::ekinToMomentum(ke_sde, m0) );
  const double mv_seb( RFlyer::ekinToMomentum(ke_seb, m0) );
  const double mv_see( RFlyer::ekinToMomentum(ke_see, m0) );
  const double mv_ext( RFlyer::ekinToMomentum(ke_ext, m0) );
  
  /// emittance steps
  const double ex_inj(0.5*RUnits::um);   double ey_inj(0.3*RUnits::um); 
  const double ex_fdb(0.5*RUnits::um);   double ey_fdb(0.3*RUnits::um);
  const double ex_fde(1.7*RUnits::um);   double ey_fde(1.1*RUnits::um);
  const double ex_feb(1.7*RUnits::um);   double ey_feb(1.1*RUnits::um); 
  const double ex_fee(0.45*RUnits::um);  double ey_fee(0.42*RUnits::um); 
  const double ex_sdb(0.45*RUnits::um);  double ey_sdb(0.42*RUnits::um); 
  const double ex_sde(2.2*RUnits::um);   double ey_sde(2.5*RUnits::um); 
  const double ex_seb(2.2*RUnits::um);   double ey_seb(2.5*RUnits::um); 
  const double ex_see(0.3*RUnits::um);   double ey_see(0.2*RUnits::um); 
  const double ex_ext(1.2*RUnits::um);   double ey_ext(0.75*RUnits::um); 
  
  /// boolean sequence checker : initialised at false = never been in sequence
  /// used to avoid computational artefact due to interpolation
  bool sc_inj(false);    /// injection
  bool sc_fdb(false);    /// first deceleration begin
  bool sc_fde(false);    /// first deceleration end
  bool sc_feb(false);    /// first electron cooling start
  bool sc_fee(false);    /// first electron cooling end
  bool sc_sdb(false);    /// second deceleration begin
  bool sc_sde(false);    /// second deceleration end
  bool sc_seb(false);    /// second electron cooling start
  bool sc_see(false);    /// second electron cooling end
  bool sc_ext(false);    /// extraction (rebunched);
  
  /// fill step array with step values
  /// time
  tcs[0] = tc_inj; tcs[1] = tc_fdb; tcs[2] = tc_fde; tcs[3] = tc_feb; 
  tcs[4] = tc_fee; tcs[5] = tc_sdb; tcs[6] = tc_sde; tcs[7] = tc_seb; 
  tcs[8] = tc_see; tcs[9] = tc_ext; 
   /// kinetic energy 
  kes[0] = ke_inj; kes[1] = ke_fdb; kes[2] = ke_fde; kes[3] = ke_feb; 
  kes[4] = ke_fee; kes[5] = ke_sdb; kes[6] = ke_sde; kes[7] = ke_seb; 
  kes[8] = ke_see; kes[9] = ke_ext; 
  /// momentum
  mvs[0] = mv_inj; mvs[1] = mv_fdb; mvs[2] = mv_fde; mvs[3] = mv_feb; 
  mvs[4] = mv_fee; mvs[5] = mv_sdb; mvs[6] = mv_sde; mvs[7] = mv_seb; 
  mvs[8] = mv_see; mvs[9] = mv_ext; 
  /// horizontal emittance
  exs[0] = ex_inj; exs[1] = ex_fdb; exs[2] = ex_fde; exs[3] = ex_feb; 
  exs[4] = ex_fee; exs[5] = ex_sdb; exs[6] = ex_sde; exs[7] = ex_seb; 
  exs[8] = ex_see; exs[9] = ex_ext; 
  /// vertical emittance
  eys[0] = ey_inj; eys[1] = ey_fdb; eys[2] = ey_fde; eys[3] = ey_feb; 
  eys[4] = ey_fee; eys[5] = ey_sdb; eys[6] = ey_sde; eys[7] = ey_seb; 
  eys[8] = ey_see; eys[9] = ey_ext; 
  /// sequence checker
  scs[0] = sc_inj; scs[1] = sc_fdb; scs[2] = sc_fde; scs[3] = sc_feb; 
  scs[4] = sc_fee; scs[5] = sc_sdb; scs[6] = sc_sde; scs[7] = sc_seb; 
  eys[8] = sc_see; scs[9] = sc_ext; 
  
  
  /// beam/machine parameters
  const double be_ncp=2e7;  /// beam - number of circulationg pbar
  
  /// sampling arrays
  double tc[ns];  /// time cycle - array of sampling point
  double ke[ns];  /// kinetic energy - array of sampling point
  double mv[ns];  /// kinetic energy - array of sampling point
  double ex[ns];  /// horizontal physical emittance (um, sigma2/TwissB)
  double ey[ns];  /// vertical physical emittance (um, sigma2/TwissB)
  
  /// time cycle
  /// ----------
  double tsv((tc_ext-tc_inj)/(ns-1));  /// time cycle - time step value
  for( int i=0; i<ns; ++i ) { tc[i] = tsv*i + tc_inj; }
  
  /// cycle : energy, emittance...
  /// ----------------------------
  /// inj : injection, initalisation
  ke[0]=kes[0];
  mv[0]=mvs[0];
  ex[0]=exs[0];
  ey[0]=eys[0];
  for( int i=1; i<ns; ++i ) { /// main loop
    for( int j=0; j<cns-1; ++j ) {  /// cycle loop
      if( tc[i]>tcs[j] && tc[i]<=tcs[j+1] ) { 
        /// mv should decrease linearly with B*rho
        mv[i] = sequ(tcs[j], tcs[j+1], tsv, mvs[j], mvs[j+1], mv[i-1], scs[j]);
        /// ke is a function of mv
        ke[i] = RFlyer::momentumToEkin(mv[i], m0);
        /// does emittance evolve like mv or ke??? linear? To confirm
        ex[i] = sequ(tcs[j], tcs[j+1], tsv, exs[j], exs[j+1], ex[i-1], scs[j]);
        ey[i] = sequ(tcs[j], tcs[j+1], tsv, eys[j], eys[j+1], ey[i-1], scs[j]);
      }
    } 
  } /// end of main loop
  
  /// reboot boolean checkers
  /*bsc_inj = false;
  bsc_fdb = false;
  bsc_fde = false;
  bsc_feb = false;
  bsc_fee = false;
  bsc_sdb = false;
  bsc_sde = false;
  bsc_seb = false;
  bsc_see = false;*/
  
  /// dump stuff
  cout << endl;
  cout << left;
  cout << setw(3)   << "i"          << "   ";
  cout << setw(10)  << "time (s)"   << "   ";
  cout << setw(10)  << "ke (keV)"   << "   ";
  cout << setw(10)  << "mv (MeV/c)" << "   ";
  cout << setw(10)  << "ex (um)"    << "   ";
  cout << setw(10)  << "ey (um)"    << "   ";
  cout << endl;
  for( int i=0; i<ns; ++i ) { 
    cout << setw(3)   << i                    << "   ";
    cout << setw(10)  << tc[i]/RUnits::s      << "   ";
    cout << setw(10)  << ke[i]/RUnits::keV    << "   ";
    cout << setw(10)  << mv[i]/RUnits::MeV_c  << "   ";
    cout << setw(10)  << ex[i]/RUnits::um     << "   ";
    cout << setw(10)  << ey[i]/RUnits::um     << "   ";
    cout << endl; 
  }
  
  /// plot stuff
  /// x pointer
  double* tcp = &tc[0];
  /// mv pointer
  double* mvp = &mv[0];
  TGraph* mvg = new TGraph(ns, tcp, mvp);
  mvg->Draw("ALP");
  gPad->Print("../results/20170705/mvg.pdf", "pdf");
  /// ke pointer
  double* kep = &ke[0];
  TGraph* keg = new TGraph(ns, tcp, kep);
  keg->Draw("ALP");
  gPad->Print("../results/20170705/keg.pdf", "pdf");
  /// emittance (h and v) pointer
  double* eyp = &ey[0];
  TGraph* eyg = new TGraph(ns, tcp, eyp);
  eyg->Draw("ALP");
  eyg->SetLineColor(2);
  double* exp = &ex[0];
  TGraph* exg = new TGraph(ns, tcp, exp);
  exg->Draw("same");
  gPad->Print("../results/20170705/emg.pdf", "pdf");
  
  
  
  /*
  /// flyer
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar, 0.037);
  //cout << endl; fly->dump();
  
  /// machine
  RMachine* mac = new RMachine("elena", 75.*RUnits::um, 3., 30.4055);
  //cout << endl; mac->dump();
  
  /// gas N2
  RGas* gas = new RGas(RParticleList::mN2(), 4.e-12*RUnits::mbar, 300.);
  cout << endl; gas->dumpRates(fly, mac);
  //delete gas;
  //gas = new RGas(RParticleList::mH2(), 4.e-12*RUnits::mbar, 300.);
  //cout << endl; gas->dumpRates(fly, mac);
  
  /// dump ionisation/cm/s
  double ion_s = npb * gas->ion_pbarnuGas(fly);
  double ion_m_s = ion_s / 30.;
  cout << endl; cout << ion_s << "  " << ion_m_s << endl;
  
  
  /// gas H2
  delete gas;
  gas = new RGas(RParticleList::mH2(), 4.e-12*RUnits::mbar, 300.);
  cout << endl; gas->dumpRates(fly, mac);
  
  /// dump ionisation/cm/s
  ion_s = npb * gas->ion_pbarnuGas(fly);
  ion_m_s = ion_s / 30.;
  cout << endl; cout << ion_s << "  " << ion_m_s << endl;
  */
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; 
  cout << "elenaCycle() function end " << endl;
  RDump::line();
}




/// elenaRate() - 22/06/2017
/// ----------------------------------------------------------------------------
/// compute elena rate with H2
void elenaRate()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;   cout << "elenaRate() - 22.06.2017" << endl;
  RDump::line();  cout << "elena rate with beam" << endl;
  
  
  /// beam parameters
  double npb=1e7;
  
  /// flyer
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar, 0.037);
  //cout << endl; fly->dump();
  
  /// machine
  RMachine* mac = new RMachine("elena", 75.*RUnits::um, 3., 30.4055);
  //cout << endl; mac->dump();
  
  /// gas N2
  RGas* gas = new RGas(RParticleList::mN2(), 4.e-12*RUnits::mbar, 300.);
  cout << endl; gas->dumpRates(fly, mac);
  //delete gas;
  //gas = new RGas(RParticleList::mH2(), 4.e-12*RUnits::mbar, 300.);
  //cout << endl; gas->dumpRates(fly, mac);
  
  /// dump ionisation/cm/s
  double ion_s = npb * gas->ion_pbarnuGas(fly);
  double ion_m_s = ion_s / 30.;
  cout << endl; cout << ion_s << "  " << ion_m_s << endl;
  
  
  /// gas H2
  delete gas;
  gas = new RGas(RParticleList::mH2(), 4.e-12*RUnits::mbar, 300.);
  cout << endl; gas->dumpRates(fly, mac);
  
  /// dump ionisation/cm/s
  ion_s = npb * gas->ion_pbarnuGas(fly);
  ion_m_s = ion_s / 30.;
  cout << endl; cout << ion_s << "  " << ion_m_s << endl;
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; 
  cout << "elenaRate() function end " << endl;
  RDump::line();
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




