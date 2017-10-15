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
#include "RBeam.h"
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












/*
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
  bsc_see = false;
  
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
  cout << "elenaCycle() function end " << endl;
  RDump::line();
}
*/










/// elenaRate() - 22/06/2017
/// ----------------------------------------------------------------------------
/// compute elena rate with H2
void elenaRate()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;		cout << "elenaRate() - 22.06.2017" << endl;
  RDump::line();	cout << "elena rate with beam" << endl;
  
  
  /// flyer
  double beta(0.037);
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar(), beta);
  //cout << endl; fly->dump();
  
  /// beam parameters
  double npb=2.e7;			/// nb of particles pre bunch
  RBeam* beam(0);
  beam = new RBeam("elenaBeam", npb, fly);
  cout << endl; cout << "dump beam " << endl; RDump::line(9, "-");
  beam->dumpLine();
  
  
  /// machine
  double A(75.*RUnits::um);				/// acceptance
  double tBeta(3.*RUnits::m);		/// twiss beta
  double len(30.4055*RUnits::m);	/// length
  RMachine* mac = new RMachine("elena", A, tBeta, len);
  //cout << endl; mac->dump();
  
  /// gas...
  double pre(4.e-12*RUnits::mbar);
  double tem(300.);
  RGas* gas(0);
  
  cout << endl;
  cout << endl; RDump::line(80, "o"); 
  cout << "molecular Nitrogen" << endl; RDump::line(80, "o"); 
  
  /// ... = pure N2
  gas = new RGas((RMolecule*)RParticleList::mN2(), pre, tem);
  cout << endl; gas->dumpRates(beam, mac); 
  /// change energy
  fly->applyEkin(100.*RUnits::keV);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac);
  /// delete gas
  delete gas; gas=0;
  
  cout << endl;
  cout << endl; RDump::line(80, "o");
  cout << "molecular hydrogen" << endl; RDump::line(80, "o"); 
  
  /// ... = pure H2
  gas = new RGas((RMolecule*)RParticleList::mH2(), pre, tem);
  
  /// change momentum
  fly->applyMomentum(100.*RUnits::MeV_c);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac); 
  /// change momentum
  fly->applyMomentum(35.*RUnits::MeV_c);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac); 
  /// change energy
  fly->applyEkin(100.*RUnits::keV);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac);
  /// delete gas
  delete gas; gas=0;
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "elenaRate() function end " << endl; RDump::line();
}










/// testRutherfordScattering() - 22.06.2017
/// ----------------------------------------------------------------------------
/// test new atom/mol/gas "old style" for rutherford scattering results
void testRutherfordScattering()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;cout << "testRutherfordScattering() - 22.06.2017" << endl; 
  RDump::line(); 
  cout << "test atom/mol/gas for rutherford scattering results" << endl;
  
  
  /// do stuff
  /// --------------------------------------------------------------------------
  /// define a gas
  double pre(4.e-12*RUnits::mbar);  /// overall pressure
  double tem(300.*RUnits::K);       /// overall temperature
  vector<RMolecule*> mov;        		/// molecules
  vector<double> mfv;               /// composition fract of molecules
  //mov.push_back((RMolecule*)RParticleList::mH2());  mfv.push_back(0.50);
  //mov.push_back((RMolecule*)RParticleList::mN2());  mfv.push_back(0.50);
  mov.push_back((RMolecule*)RParticleList::mN2());  mfv.push_back(1.00);
  RGas* gas = new RGas("elena test gas", pre, tem);
  gas->setCompound(mov, mfv);
  /// dump the gas
  cout << endl; gas->dump(cout, "pkr261");
  
  /// define machine
  double len(30.4055*RUnits::m);	/// machine length
  double tbp(3.*RUnits::m);				/// twiss beta parameter at detector
  double tap(75.*RUnits::um);			/// twiss acceptance param = tbpMax/tbp
  RMachine* mac(0);
  mac = new RMachine("elena", tap, tbp, len);
  /// dump the machine
  cout << endl; 
  cout << endl; cout << "dump machine parameters" << endl; RDump::line(23, "-");
  mac->dumpLine();

  /// define flyers
  double beta(0.037);
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar(), beta);
  /// dump the flyer
  cout << endl; 
  cout << endl; cout << "dump flyer parameters" << endl; RDump::line(21, "-");
  fly->dumpLine(); fly->dumpBeta();
  
  /// atom dumpcs test
  RAtom* ato = gas->getCompound()[0]->getComposition()[0];
  cout << endl; 
  cout << endl; ato->dumpcs(fly, mac);
  
  /// molecule dumpcs test
  cout << endl; 
  for( int i=0; i<gas->getCompound().size(); ++i ) {
    cout << endl; gas->getCompound()[i]->dumpcs(fly, mac);
  }
  
  /// dump rates
  cout << endl; 
  cout << endl; gas->dumpRates(fly, mac);
  delete gas; gas = new RGas((RMolecule*)RParticleList::mH2(), pre, tem);
  cout << endl; 
  cout << endl; gas->dumpRates(fly, mac);
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  delete gas; gas=0; delete fly; fly=0; delete mac; mac=0;
  cout << endl; cout << "testRutherfordScattering() end " << endl; 
  RDump::line();
}









/// testMacFly() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Machine/Flyer classs
void testMacFly()
{
  /// function header
  cout << endl; cout << "testMacFly() - 15.06.2017" << endl; RDump::line(); 
  cout << "test Machine/Flyer class" << endl;
  
  
  /// middle
  RParticle* par(0);
  RFlyer* fly(0);
  double beta(0);
  double kin(0);
  double mv(0);
  par = (RParticle*)RParticleList::hpbar();
  fly = new RFlyer(par, beta); delete par; par=0;
  
  cout << "the flyer is: ";
  fly->dumpLine();
  cout << endl; cout << "beta=" << beta << " gives: "<< endl;
  fly->dumpBeta();
  
 	kin=5.3*RUnits::MeV; fly->applyEkin(kin);
  cout << endl; cout << "ekin=" << kin/RUnits::MeV << "MeV" << " gives: ";
  cout << endl; fly->dumpBeta();
  
  mv=35.*RUnits::MeV_c; fly->applyMomentum(mv);
 	cout << endl; cout << "mv=" << mv/RUnits::MeV_c << "MeV/c" << " gives: ";
 	cout << endl; fly->dumpBeta();
  
  kin=100.*RUnits::keV; fly->applyEkin(kin);
  cout << endl; cout << "ekin=" << kin/RUnits::MeV << "MeV" << " gives: ";
  cout << endl; fly->dumpBeta();
 
  
  /// end of function
  delete fly; fly=0;
  cout << endl; cout << "tMac() function end " << endl; RDump::line(); 
}





/// testGas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class, especially dump and dumpRates
void testGas()
{
  /// function header
  cout << endl; cout << "testGas() - 15.06.2017" << endl; RDump::line(); 
  cout << "test Gas class function start" << endl;
  
  
  /// gas properties
  double pre(4.e-12*RUnits::mbar);  /// overall pressure
  double tem(298.*RUnits::K);       /// overall temperature
  RGas* gas(0);											/// gas pointer
  /// define a gas
  vector<RMolecule*> mov(0);        /// molecules vector
  vector<double> mfv;               /// composition fract of molecules
  mov.push_back((RMolecule*)RParticleList::mH2());  mfv.push_back(0.91);
  mov.push_back((RMolecule*)RParticleList::mCO2()); mfv.push_back(0.05);
  mov.push_back((RMolecule*)RParticleList::mCO());  mfv.push_back(0.03);
  mov.push_back((RMolecule*)RParticleList::mCH4()); mfv.push_back(0.01);
  gas = new RGas("elenaVac", pre, tem);
  gas->setCompound(mov, mfv);
  
  /// dump the gas
	cout << endl << "test elena gas with molecules mixing" << endl;
  RDump::line(36, "-");
  cout << endl; gas->dump(cout, "pkr261"); delete gas; gas = 0;
  
  /// simple H2 gas
  gas = new RGas((RMolecule*)RParticleList::mH2(), pre, tem);
  cout << endl << "test simple H2 gas " << endl;
  RDump::line(18, "-");
  cout << endl; gas->dump(cout, "pkr261"); delete gas; gas = 0;
  
  
  /// end of function
  cout << endl; cout << "testGas() end" << endl;  RDump::line(); 
}









/// testval_csion_pbarGas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Ionisation cross section eval/interpolation function
void testval_csion_pbarGas()
{
  /// function header
  cout << endl; cout << "testval_csion_pbarGas() - 15.06.2017" << endl; 
  RDump::line(); cout << "check pbar ionisation cs evaluation" << endl;
  
  /// do stuff 
  cout << (RParticleList::mCO2())->ion_pbarcs(500.*RUnits::keV) 
       << "m2" << endl;
  
  /// end of function
  cout << endl; cout << "testval_csion_pbarGas()" << endl; RDump::line(); 
}







/// plotIon_pbarcs() - 14.06.2017
/// ----------------------------------------------------------------------------
/// plot ion_pbarcs data graphs
void plotIon_pbarcs(string const& today)
{
  /// function header
  cout << endl; cout << "plotIon_pbarcs(string today) - 14.06.2017" << endl; 
  RDump::line(); cout << "plot ion_pbarcs data graphs" << endl;
  
  
  /// do stuff
  /// string manip
  string pPrefix("../results/");
  string path(pPrefix+today+"/");
  string fPrefix("cs_ion_pbarData_");
  string file("");
  string fSuffix(".pdf");
  /// molecule and graph
  const RMolecule* mol(0);
  TGraphErrors* tge(0);
  
  /// molecular argon
  mol = RParticleList::mAr(); tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mAr.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular ch4
  mol = RParticleList::mCH4();  tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mCH4.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular CO
  mol = RParticleList::mCO();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  delete tge; delete mol; delete gPad;
  
  /// molecular CO2
  mol = RParticleList::mCO2();  tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mCO2.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular H2
  mol = RParticleList::mH2();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mH2.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular krypton
  mol = RParticleList::mKr();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mKr.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular nitrogen
  mol = RParticleList::mN2();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mN2.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular neon
  mol = RParticleList::mNe();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mNe.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular xenon
  mol = RParticleList::mXe();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mXe.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  
  /// end of function
  tge = 0; mol = 0; gPad = 0;
  cout << endl; cout << "plotIon_pbarcs() end" << endl; RDump::line();
}






/// testMolecule() - 14.06.2017
/// ----------------------------------------------------------------------------
/// test RMolecule class
void testMolecule()
{
  /// function header
  cout << endl << "testMolecule() - 14.06.2017" << endl; RDump::line();
  cout << "test Molecule class function start" << endl;
  
  
  /// do stuff
  RMolecule* mol(0);
  
  /// syntax 1
	vector<RAtom*> mcv;	/// molecule composition vector
  mcv.push_back((RAtom*)RParticleList::aC12());
  mcv.push_back((RAtom*)RParticleList::aO16());
  mcv.push_back((RAtom*)RParticleList::aO16());
  mol = new RMolecule("Carbon dioxide", "CO2", "C12-O16_2", mcv);
	cout << endl; cout << "syntax 1:" << endl;
  mol->dump(); delete mol; mol = 0;
  
  /// syntax 2
  mol = new RMolecule("repMol", "rM", "r30M0");
  mol->addAtom((RAtom*)RParticleList::aN14());
  mol->addAtom((RAtom*)RParticleList::aH1());
  mol->addAtom((RAtom*)RParticleList::aH1());
  mol->addAtom((RAtom*)RParticleList::aH1());
  cout << endl; cout << "syntax 2:" << endl;
  mol->dump(); delete mol; mol = 0;
  
  /// syntax 3
  mol = (RMolecule*)RParticleList::mCH4();
 	cout << endl; cout << "syntax 3:" << endl;
  mol->dump(); delete mol; mol = 0;
  
  
  /// end of function
  cout << endl; cout << "testMolecule() end" << endl; RDump::line();
}





/// listParticles() - 07.07.2017
/// ----------------------------------------------------------------------------
/// list all particles in particleList namespace
void listParticles()
{
  /// function header
  cout << endl << "listParticles() - 07.07.2017" << endl; RDump::line();
  cout << "list particles in the namespace RParticleList.cc" << endl;
  
  const RAtom* atom(0);
  cout << endl; 
  cout << "dump atoms (line)" << endl;
  RDump::line(21, "-");
  atom = RParticleList::aH1();     atom->dumpLine(); delete atom; 
  atom = RParticleList::aC12();    atom->dumpLine(); delete atom;
  atom = RParticleList::aN14();    atom->dumpLine(); delete atom;
  atom = RParticleList::aO16();    atom->dumpLine(); delete atom;
  atom = RParticleList::aNe20();   atom->dumpLine(); delete atom;
  atom = RParticleList::aAr40();   atom->dumpLine(); delete atom;
  atom = RParticleList::aKr84();   atom->dumpLine(); delete atom;
  atom = RParticleList::aXe132();  atom->dumpLine(); delete atom;
  atom = 0;
  
  const RMolecule* mol(0);
  cout << endl << endl;
  cout << "dump molecules" << endl;
  RDump::line(14, "-");
  mol = RParticleList::mAr();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCH4();  mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCO();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCO2();  mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mH2();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mKr();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mN2();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mNe();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mXe();   mol->dumpLine(); cout << endl; delete mol;
  mol = 0;
  
  const RParticle* par(0);
  cout << endl << endl;
  cout << "dump hadrons" << endl;
  RDump::line(12, "-");
  par = RParticleList::hpbar();   par->dumpLine(); cout << endl; delete par;
  par = 0;
  
  
  /// end of function
  cout << endl; cout << "listParticles() function end" << endl; RDump::line();
}






/// testAtom() - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RAtom class
void testAtom()
{
  /// function header
  cout << "testAtom() - 13.06.2017" << endl; RDump::line();
  cout << "test Atom class function start" << endl;
  
  
  /// test class, basics
  RAtom* ato1 = new RAtom("repAtom", "Re", 100., 200.);
  cout << endl; ato1->dump();
  
  cout << endl << "test atom get accessors" << endl;
  RDump::line(23, "-");
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

/// testParticle() - 12.06.2017
/// ----------------------------------------------------------------------------
/// test class RParticle
void testParticle()
{
  /// function header
  cout << "testParticle() - 12.06.2017" << endl; RDump::line();
  cout << "test RParticle class function start" << endl;
  
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
  cout << endl; cout << "testParticle() function end" << endl; RDump::line();
}




/// testObject() - 09.06.2017
/// ----------------------------------------------------------------------------
/// test class RObject
void testObject()
{
  /// function header
  cout << "testObject() - 09.06.2017" << endl; RDump::line();
  cout << "test RObject class function start" << endl;
  
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
  
  /// end of function
  cout << endl; cout << "testObject() function end" << endl; RDump::line();
}




