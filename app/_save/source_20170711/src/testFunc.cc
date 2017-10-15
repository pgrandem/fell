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
#include <fstream>
//#include <iomanip>
/// root classes
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TObject.h"
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
#include "RMath.h"
#include "RParticleList.h"
#include "RUnits.h"
/// rep data functions
#include "elenaCool.h"
/// local functions
#include "testFunc.h"
using namespace std;
//using namespace RParticleList;
//using namespace RUnits;

/// category 01
/// ****************************************************************************

/// ecplots()
/// ----------------------------------------------------------------------------
/// elena cycle plots - mv, ke, ex and ey
void ecplots()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;   cout << "ecplots() - 06/07/2017" << endl;
  RDump::line();  cout << "ionistations/m/s as function of time" << endl;
  
  
  /// do the stuff
  /// --------------------------------------------------------------------------
  int const np(53);         /// number of points
  TGraph* mvg = ecmvg(np);   /// momentum cycle
  TGraph* keg = eckeg(np);   /// energy cycle
  TGraph* exg = ecxeg(np);   /// x emittance cycle
  TGraph* eyg = ecyeg(np);   /// y emittance cycle
  
  /// dump stuff
  double* mvx = mvg->GetX();
  double* mvy = mvg->GetY();
  double* key = keg->GetY();
  double* exy = exg->GetY();
  double* eyy = eyg->GetY();
  cout << endl;
  cout << left;
  cout << setw(3)   << "i"          << "   ";
  cout << setw(10)  << "time (s)"   << "   ";
  cout << setw(10)  << "mv (MeV/c)" << "   ";
  cout << setw(10)  << "ke (MeV)"   << "   ";
  cout << setw(10)  << "ex (um)"    << "   ";
  cout << setw(10)  << "ey (um)"    << "   ";
  cout << endl;
  for( int i=0; i<np; ++i ) { 
    cout << setw(3)   << i      << "   ";
    cout << setw(10)  << mvx[i] << "   ";
    cout << setw(10)  << mvy[i] << "   ";
    cout << setw(10)  << key[i] << "   ";
    cout << setw(10)  << exy[i] << "   ";
    cout << setw(10)  << eyy[i] << "   ";
    cout << endl; 
  }
  
  /// draw stuff
  mvg->SetMarkerStyle(4);
  mvg->SetMarkerSize(0.4);
  mvg->Draw("ALP");
  gPad->Print("../results/20170706/mvg.pdf", "pdf");
  keg->SetMarkerStyle(4);
  keg->SetMarkerSize(0.4);
  keg->Draw("ALP");
  gPad->Print("../results/20170706/keg.pdf", "pdf");
  exg->SetMarkerStyle(4);
  exg->SetMarkerSize(0.4);
  exg->Draw("ALP");
  gPad->Print("../results/20170706/exg.pdf", "pdf");
  eyg->SetMarkerStyle(4);
  eyg->SetMarkerSize(0.4);
  eyg->Draw("ALP");
  gPad->Print("../results/20170706/eyg.pdf", "pdf");
  
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "ecplots() function end " << endl; RDump::line();
}

/// ecxeg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle x emittance graph (x=horizontal)
TGraph* ecxeg(int const np)
{ 
  /// cycle steps 
  const int cns(10);  /// cycle number of steps
  double xsa[cns];    /// x steps array
  double ysa[cns];    /// y steps array
  
  /// time steps - tc = time cycle
  const double xUnit = RUnits::s;
  const double tc_inj( 0.*xUnit); /// injection;
  const double tc_fdb( 1.*xUnit); /// first deceleration begin
  const double tc_fde( 5.*xUnit); /// first deceleration end
  const double tc_feb( 6.*xUnit); /// first electron cooling start
  const double tc_fee(16.*xUnit); /// first electron cooling end
  const double tc_sdb(17.*xUnit); /// second deceleration begin
  const double tc_sde(21.*xUnit); /// second deceleration end
  const double tc_seb(22.*xUnit); /// second electron cooling start
  const double tc_see(25.*xUnit); /// second electron cooling end
  const double tc_ext(26.*xUnit); /// extraction (rebunched);
  
  /// x emittance steps
  const double yUnit = RUnits::um;
  const double inj(0.5*yUnit);
  const double fdb(0.5*yUnit);
  const double fde(1.7*yUnit);
  const double feb(1.7*yUnit);
  const double fee(0.45*yUnit);
  const double sdb(0.45*yUnit);
  const double sde(2.2*yUnit);
  const double seb(2.2*yUnit);
  const double see(0.3*yUnit);
  const double ext(1.2*yUnit);
  
  /// attribute values
  xsa[0]=tc_inj;  xsa[1]=tc_fdb;  xsa[2]=tc_fde;  xsa[3]=tc_feb; 
  xsa[4]=tc_fee;  xsa[5]=tc_sdb;  xsa[6]=tc_sde;  xsa[7]=tc_seb; 
  xsa[8]=tc_see;  xsa[9]=tc_ext; 
  
  ysa[0]=inj;   ysa[1]=fdb;   ysa[2]=fde;   ysa[3]=feb; 
  ysa[4]=fee;   ysa[5]=sdb;   ysa[6]=sde;   ysa[7]=seb; 
  ysa[8]=see;   ysa[9]=ext;
  
  /// output
  TGraph* tg = cg(np, cns, xsa, ysa, "cooling");
  double* x = tg->GetX();
  double* y = tg->GetY();
  for( int i=0; i<tg->GetN(); ++i ) { 
    tg->SetPoint(i, x[i]/xUnit, y[i]/yUnit);
  }
  tg->SetName("ex2t");
  tg->SetTitle("elena cycle - horizontal emittance");
  tg->GetXaxis()->SetTitle("time  (s)");
  tg->GetYaxis()->SetTitle("horizontal emittance  (um)");
  return tg;
}

/// ecyeg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle y emittance graph (x=horizontal)
TGraph* ecyeg(int const np)
{ 
  /// cycle steps 
  const int cns(10);  /// cycle number of steps
  double xsa[cns];    /// x steps array
  double ysa[cns];    /// y steps array
  
  /// time steps - tc = time cycle
  const double xUnit = RUnits::s;
  const double tc_inj( 0.*xUnit); /// injection;
  const double tc_fdb( 1.*xUnit); /// first deceleration begin
  const double tc_fde( 5.*xUnit); /// first deceleration end
  const double tc_feb( 6.*xUnit); /// first electron cooling start
  const double tc_fee(16.*xUnit); /// first electron cooling end
  const double tc_sdb(17.*xUnit); /// second deceleration begin
  const double tc_sde(21.*xUnit); /// second deceleration end
  const double tc_seb(22.*xUnit); /// second electron cooling start
  const double tc_see(25.*xUnit); /// second electron cooling end
  const double tc_ext(26.*xUnit); /// extraction (rebunched);
  
  /// y emittance steps
  const double yUnit = RUnits::um;
  const double inj(0.3*yUnit);
  const double fdb(0.3*yUnit);
  const double fde(1.1*yUnit);
  const double feb(1.1*yUnit);
  const double fee(0.42*yUnit);
  const double sdb(0.42*yUnit);
  const double sde(2.5*yUnit);
  const double seb(2.5*yUnit);
  const double see(0.2*yUnit);
  const double ext(0.75*yUnit);
  
  /// attribute values
  xsa[0]=tc_inj;  xsa[1]=tc_fdb;  xsa[2]=tc_fde;  xsa[3]=tc_feb; 
  xsa[4]=tc_fee;  xsa[5]=tc_sdb;  xsa[6]=tc_sde;  xsa[7]=tc_seb; 
  xsa[8]=tc_see;  xsa[9]=tc_ext; 
  
  ysa[0]=inj;   ysa[1]=fdb;   ysa[2]=fde;   ysa[3]=feb; 
  ysa[4]=fee;   ysa[5]=sdb;   ysa[6]=sde;   ysa[7]=seb; 
  ysa[8]=see;   ysa[9]=ext;
  
  /// output
  TGraph* tg = cg(np, cns, xsa, ysa, "cooling");
  double* x = tg->GetX();
  double* y = tg->GetY();
  for( int i=0; i<tg->GetN(); ++i ) { 
    tg->SetPoint(i, x[i]/xUnit, y[i]/yUnit);
  }
  tg->SetName("ex2t");
  tg->SetTitle("elena cycle - vertical emittance");
  tg->GetXaxis()->SetTitle("time  (s)");
  tg->GetYaxis()->SetTitle("vertical emittance  (um)");
  return tg;
}



/// ecekg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle kinetic energy graph 
/// carefull: ecmvg has a unit system "hardcoded"
/// -> assuming emittance growth is proportional to pomentum
TGraph* eckeg(int const np)
{
  TGraph* keg = new TGraph(np);
  TGraph* mvg = ecmvg(np);
  double* mvx = mvg->GetX();
  double* mvy = mvg->GetY();
  const double m0 = 1.0072765*RUnits::amu;
  for( int i=0; i<np; ++i ) {
    /// warning : unit system!
    double key( RFlyer::momentumToEkin(mvy[i]*RUnits::MeV_c, m0)/RUnits::MeV );
    keg->SetPoint(i, mvx[i] , key);
  }
  /// output
  keg->SetName("ke2t");
  keg->SetTitle("elena cycle - kinetic energy");
  keg->GetXaxis()->SetTitle("time  (s)");
  keg->GetYaxis()->SetTitle("kinetic energy  (MeV)");
  return keg;
}



/// ecmvg()
/// ----------------------------------------------------------------------------
/// elena cycle momentum graph
TGraph* ecmvg(int const np)
{ 
  /// cycle steps 
  const int cns(10);  /// cycle number of steps
  double xsa[cns];    /// x steps array
  double ysa[cns];    /// y steps array
  
  /// time steps - tc = time cycle
  const double xUnit = RUnits::s;
  const double tc_inj( 0.*xUnit); /// injection;
  const double tc_fdb( 1.*xUnit); /// first deceleration begin
  const double tc_fde( 5.*xUnit); /// first deceleration end
  const double tc_feb( 6.*xUnit); /// first electron cooling start
  const double tc_fee(16.*xUnit); /// first electron cooling end
  const double tc_sdb(17.*xUnit); /// second deceleration begin
  const double tc_sde(21.*xUnit); /// second deceleration end
  const double tc_seb(22.*xUnit); /// second electron cooling start
  const double tc_see(25.*xUnit); /// second electron cooling end
  const double tc_ext(26.*xUnit); /// extraction (rebunched);
  
  /// momentum steps - mv = momentum (mass*velocity)
  /// should correspond to magnetic cycle steps:
  /// -> B*rho # momentum (# = proportional)
  /// shoud be a linear dependance of time
  const double yUnit = RUnits::MeV_c;
  const double mv_inj( 100.*yUnit );  
  const double mv_fdb( 100.*yUnit );  
  const double mv_fde(  35.*yUnit );
  const double mv_feb(  35.*yUnit );  
  const double mv_fee(  35.*yUnit );
  const double mv_sdb(  35.*yUnit );
  const double mv_sde( 13.7*yUnit );
  const double mv_seb( 13.7*yUnit );
  const double mv_see( 13.7*yUnit );
  const double mv_ext( 13.7*yUnit );
  /// attribute values
  xsa[0]=tc_inj;  xsa[1]=tc_fdb;  xsa[2]=tc_fde;  xsa[3]=tc_feb; 
  xsa[4]=tc_fee;  xsa[5]=tc_sdb;  xsa[6]=tc_sde;  xsa[7]=tc_seb; 
  xsa[8]=tc_see;  xsa[9]=tc_ext; 
  
  ysa[0] = mv_inj; ysa[1] = mv_fdb; ysa[2] = mv_fde; ysa[3] = mv_feb; 
  ysa[4] = mv_fee; ysa[5] = mv_sdb; ysa[6] = mv_sde; ysa[7] = mv_seb; 
  ysa[8] = mv_see; ysa[9] = mv_ext;
  
  /// output
  TGraph* tg = cg(np, cns, xsa, ysa, "linear");
  double* x = tg->GetX();
  double* y = tg->GetY();
  for( int i=0; i<tg->GetN(); ++i ) { 
    tg->SetPoint(i, x[i]/xUnit, y[i]/yUnit);
  }
  tg->SetName("p2t");
  tg->SetTitle("elena cycle - momentum");
  tg->GetXaxis()->SetTitle("time  (s)");
  tg->GetYaxis()->SetTitle("momentum  (MeV/c)");
  return tg;
}

/// cg(int ng, int ns, double xs[], double ys[], std::string ymodel);
/// ----------------------------------------------------------------------------
/// cycle graph from array of steps and model to interpol y(x) in between steps
/// ng: nb of graph points
/// ns: nb of steps points
/// xs: x steps array
/// ys: y steps array
/// ymodel: model to interpolate y(x) when xi < x < xi+1
TGraph* cg(int np, int cns, double xsa[], double ysa[], string const& ymodel)
{
  /// sampling arrays
  double gxa[np];  /// graph x array
  double gya[np];  /// graph y array
  double xstep( (xsa[cns-1]-xsa[0])/(np-1) );
  
  /// init first point
  gxa[0] = xsa[0];
  gya[0] = ysa[0];
  
  
  /// main loops for the different model of interpolation
  /// ---------------------------------------------------
  
  /// linear model : linear interpolation
  if( ymodel=="linear" ) {   
    for( int i=1; i<np; ++i ) { /// loop over sampling points
      /// x array
      gxa[i] = xstep*i + xsa[0];
      /// y array
      for( int j=0; j<cns-1; ++j ) {  /// loop over cycle step points
        if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] ) { 
          gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
      }
    } 
  }/// end of linear case
    
  /// cooling model : linear at deceleration, cooling curve at cooling
  else if( ymodel=="cooling" ) {
    /// scaled cooling graph (betaCool sim)
    int const fcs(3); /// first cooling step
    int const scs(7); /// second cooling step
    TGraph* fcg = scool("35MeV_c", xsa[fcs], ysa[fcs], ysa[fcs+1]); 
    TGraph* scg = scool("100keV",  xsa[scs], ysa[scs], ysa[scs+1]); 
    for( int i=1; i<np; ++i ) { /// loop over sampling points
      /// x array
      gxa[i] = xstep*i + xsa[0];
      /// y array
      for( int j=0; j<cns-1; ++j ) {  /// loop over cycle step points
        if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] && j==fcs ) { /// 1st cooling
          gya[i]=vcool(fcg, gxa[i]);
          //gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
        else if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] && j==scs ) { /// second cooling
          gya[i]=vcool(scg, gxa[i]);
          //gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
        else if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] ) { /// linear
          gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
      }
    } 
    delete fcg; fcg=0;
    delete scg; scg=0;
  } /// end of cooling case
  
  /// output!
  TGraph* tg = new TGraph(np);
  for( int i=0; i<np; ++i ) { 
    tg->SetPoint(i, gxa[i], gya[i]);
  }
  /// end of routine
  return tg;
}







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
  if( bsc==false ) { init = vsi; bsc=true; /*cout << bsc << endl;*/ }
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
  const int ns(51);     /// number of samples
  
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
  mv[0]=mvs[0];
  ke[0]=kes[0];
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



/*
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
*/


/*
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
*/






/*
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
  *//*
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "tInt() function end " << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl;
}
*/














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




/*
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
*/







/*
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

*/






/*
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
*/

/*

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
*/



/*
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
  *//*
  
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
  *//*
  
  
  
  /// end of function
  cout << endl; cout << "test Molecule class \"tMol\" function end" << endl;
  for( int i=0; i<80; ++i) { cout << "*"; } cout << endl; /// 80 "-" line 
}

*/

/// listParticles() - 07.07.2017
/// ----------------------------------------------------------------------------
/// list all particles in particleList namespace
void listParticles()
{
  /// function header
  cout << endl << "tPList()" << endl;
  RDump::line();
  cout << "test ParticleList namespace with Atom object declaration" << endl;
  
  const RAtom* atom(0);
  cout << endl; 
  cout << "dump all atoms (line)" << endl;
  RDump::line(20, "-");
  atom = RParticleList::aH1();     atom->dumpLine(); delete atom; 
  atom = RParticleList::aC12();    atom->dumpLine(); delete atom;
  atom = RParticleList::aN14();    atom->dumpLine(); delete atom;
  atom = RParticleList::aO16();    atom->dumpLine(); delete atom;
  atom = RParticleList::aNe20();   atom->dumpLine(); delete atom;
  atom = RParticleList::aAr40();   atom->dumpLine(); delete atom;
  atom = RParticleList::aKr84();   atom->dumpLine(); delete atom;
  atom = RParticleList::aXe132();  atom->dumpLine(); delete atom;
  
  const RMolecule* mol(0);
  cout << endl;
  cout << "dump all molecules" << endl;
  RDump::line(60, "-");
  mol = RParticleList::mAr();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCH4();  mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCO();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCO2();  mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mH2();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mKr();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mN2();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mNe();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mXe();   mol->dumpLine(); cout << endl; delete mol;
  
  
  
  /// end of function
  cout << endl; cout << "tAtom() function end" << endl; RDump::line();
}




/// testAtom - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RAtom class
void tAtom()
{
  /// function header
  cout << "tAtom()" << endl;
  RDump::line(); /// 80 "-" line 
  cout << "test Atom class function start" << endl;
  
  
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
  cout << endl; cout << "tAtom() function end" << endl; RDump::line();
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
  RParticle* par1(0); 
  RParticle* par2(0); 
  RParticle* par3(0); 
  par1 = new RParticle("proton", "p", 1.*RUnits::eV, 938.*RUnits::MeV_c2); 
  par2 = new RParticle("neutron", "n", 0., 940.*RUnits::MeV_c2); 
  par3 = new RParticle("Helium 2+", "He2+", 2.*RUnits::eV, 4.*RUnits::MeV_c2); 
  
  cout << endl; testPol(obj1);
  cout << endl; testPol(par1);
  cout << endl; testPol(par2);
  cout << endl; testPol(par3);
  
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
  
  /// end of function
  delete par1; par1=0;
  delete par2; par2=0;
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




