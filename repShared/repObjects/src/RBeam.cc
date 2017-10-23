/// RBeam.cc
/// ----------------------------------------------------------------------------
/// rep RAtom c++ class 
/// Pierre GRANDEMANGE
/// 13/07/2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <string>
/// root classes
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"
/// rep namespaces
#include "RDump.h"
#include "RMath.h"
#include "RParticleList.h"
#include "RUnits.h"
/// rep classes
#include "RBeam.h"
#include "RFlyer.h"
#include "RGas.h"
#include "RMachine.h"
using namespace std;


/// constructor, destructor, copy 
/// ****************************************************************************

RBeam::RBeam(string const& name, int npb, RParticle* par, string const dtype[6],
             double mv0, double emit[3], double macPar[3]) : 
  RObject(name), 
  rflyNumb(npb), 
  rparType(par),
	rtitl{ "x", "y", "z", "px", "py", "pz" },
  rtype{ dtype[0], dtype[1], dtype[2], dtype[3], dtype[4], dtype[5] },
  remit{ emit[0], emit[1], emit[2] },
  runit{ RUnits::mm,    RUnits::mm,    RUnits::mm, 
         RUnits::MeV_c, RUnits::MeV_c, RUnits::MeV_c }
{
  /// allocate pointers
  rflyColl = new RFlyer[rflyNumb];
  /// build the beam
  sizeFromEmit(mv0, macPar); /// compute size from ref momentum and emittance
  func();     /// create the normalised distribution function
  flyColl();  /// random flyer collection along distribution
  hist(100);  /// fill histogram with flyer collection
  funcNH();   /// norm dist func to hist.integral to superimpose
}

RBeam::RBeam(string const& name, int npb, RParticle* par, 
             string const dtype[6], double mean[6], double size[6]) : 
	RObject(name), 
  rflyNumb(npb), 
  rparType(par),
  rmean{ mean[0], mean[1], mean[2], mean[3], mean[4], mean[5] },
  rtitl{ "x", "y", "z", "px", "py", "pz" },
  rtype{ dtype[0], dtype[1], dtype[2], dtype[3], dtype[4], dtype[5] },
  rsize{ size[0], size[1], size[2], size[3], size[4], size[5] },
  remit{  1., 1., 1. }/*compute emittance ?!!*/,
  runit{ RUnits::mm,     RUnits::mm,    RUnits::mm, 
         RUnits::MeV_c,  RUnits::MeV_c, RUnits::MeV_c }
{
  /// allocate pointer to array of flyers...
  rflyColl = new RFlyer[rflyNumb];
  /// build the beam
  func();     /// create the normalised distribution function
  flyColl();  /// build random flyColl
  hist(100);  /// fill histogram with flyer collection
  funcNH();   /// norm dist func to hist.integral to superimpose
}

RBeam::~RBeam()
{ 
  if(rflyColl!=NULL) { delete[] rflyColl; rflyColl=0; } 
  delete rparType; rparType=0;
}





/// accessors
/// ****************************************************************************





/// dump methods
/// ****************************************************************************
/// dumpFlyColl
/// --------------------------------------------------------------------------
/// dump flyer collection x y z px py pz 
void RBeam::dumpFlyColl(int modulo, ostream& flux) const
{
  /// save previous stream format settings
  ios::fmtflags f(flux.flags());        /// save current flux flags
  streamsize    p(flux.precision());    /// save current flux precision
  /// new stream format... to dump iterator i value
  flux.flags(ios::left | ios::scientific); /// ajust on right, scientifi
  const string sep(" ");
  const int iw(5);
  flux.precision(0);
  /// dump header
  RDump::header("RBeam::dumpFlyColl");
  /// dump stuff
  flux << setw(iw) << "i" << sep;
  flux << setw(10) << " x" << sep;
  flux << setw(10) << " y" << sep;
  flux << setw(10) << " z" << sep;
  flux << setw(13) << " px" << sep;
  flux << setw(13) << " py" << sep;
  flux << setw(13) << " pz";
  flux << "\n";
  //RDump::line(78, "-");
  for( int i=0; i<rflyNumb; ++i ) {
    if( i % modulo == 0 ) {
      flux << setw(iw) << (double)i << sep;
      rflyColl[i].line6D(flux);
    }
  }
  /// restore previous stream format settings
  flux.flags(f);        /// restore flags
  flux.precision(p);    /// restore precision
}



//// dumpIntegrals
/// ----------------------------------------------------------------------------
/// dump different integral values
void RBeam::dumpIntegrals(std::ostream& flux) const
{
  RDump::header("RBeam::dumpIntegrals");
  for( int i=0; i<6; ++i ) {
    TF1   f( this->getFunc()[i] );
    TH1D  h( this->getHist()[i] );
    double xM( f.GetXmax() );
    double xm( f.GetXmin() );
    cout << "h[i].Int()="  << left << setw(6) << h.Integral()        << " ";
    cout << "h[i].Int(w)=" << left << setw(7) << h.Integral("width") << " ";
    cout << "f[i].Int()="  << left << setw(8) << f.Integral(xm, xM)  << endl;
  }
}



/// dumpLine
/// --------------------------------------------------------------------------
/// dump main beam properties on one line
void RBeam::dumpLine(ostream &flux) const
{
	
	ios::fmtflags f(flux.flags());	/// save current flags in flux
  
  flux << this->getName() << " | ";
  flux << "npb=" << this->getN() << " | ";
  //flux.precision(2); flux << scientific;
  //flux << "flyer: " << this->getFlyer()->getSymbol() << " | ";
  //flux << "beta=" << this->getFlyer()->beta() << " | ";
  flux << endl; flux.flags(f);		/// restore flags
}



// dumpSize
/// ----------------------------------------------------------------------------
/// dump mean +- size (1sigma) | disTyp for x y z px py pz
void RBeam::dumpSize(std::ostream& flux) const
{
  /// save flux previous format settings
  ios::fmtflags pf(flux.flags());        /// save flux previous flags
  streamsize    pp(flux.precision());    /// save flux previous precision
  
  /// new current flux format
  const std::string sep("  ");       /// separator
  const int wid(10);
  /// dump header
  RDump::header("RBeam::dumpSize");
  /// loop over flyers position coordinate
  for( int i=0; i<6; ++i ) {
    flux << left << setw(2) << rtitl[i] << sep << rtype[i] << sep;
    flux << right << rmean[i]/runit[i] << "+-" << rsize[i]/runit[i];
    flux << runit[i].s() << '\n';
  }
  
  /// restore flux previous format settings
  flux.flags(pf);        /// restore flux previous flags
  flux.precision(pp);    /// restore flux previous precision
}
  

/// pdf6D
/// ----------------------------------------------------------------------------
/// save the beam 6 coordinates distribution function/histogram in pdf
/// returns pdf path
std::string RBeam::pdf6D() const
{
  RDump::header("RBeam::pdf6D");
  TCanvas* can(RDump::canzenbook3("can", "pdf6D"));
  const std::string pdf(RDump::pathNow(this->getName()+"_dist.pdf"));
  const std::string pdo(pdf+"[");
  const std::string pdc(pdf+"]");
  can->Print(pdo.c_str());
  for( int i=0; i<6; ++i ) {
    this->getFunc()[i].Draw(); 
    this->getHist()[i].Draw("same"); 
    can->Print(pdf.c_str()); 
  }
  can->Print(pdc.c_str());
  delete can; can=0;
  return pdf;
}




/// other methods
/// ****************************************************************************

/// flyColl
/// ----------------------------------------------------------------------------
/// Set rflyColl from rfunc and rflyNumb
/// Look at RMath::random
void RBeam::flyColl(int ranSeed)
{ 
  double *p[6];  /// array of {double collection pointer}
  for( int i=0; i<6; ++i ) {  /// loop to compute double collec.
    p[i] = RMath::reprand( &rfunc[i], runit[i], rflyNumb, ranSeed );
  }
  for( int i=0; i<rflyNumb; ++i ) { /// dbl loop to assign to flyers.
    //double d[6]={ p[0][i]*runit[0], p[1][i]*runit[1], p[2][i]*runit[2], 
    //              p[3][i]*runit[3], p[4][i]*runit[4], p[5][i]*runit[5] };
    double d[6]={ p[0][i], p[1][i], p[2][i], p[3][i], p[4][i], p[5][i] };
    rflyColl[i] = RFlyer(rparType, d);
    /// verbose/debug
    //rflyColl[i].line6D();
  }
  /// release memory
  for( int i=0; i<6; ++i ) { delete[] p[i]; p[i]=0; }
}

/// flyColl1D
/// ----------------------------------------------------------------------------
/// returns a double ptr to 1D(over 6) of a flyer collection
///   - whichD | which dim to return (0 to 5)
double*  RBeam::flyColl1D( int whichD) const
{
  double *out = new double[rflyNumb];
  for( int i=0; i<rflyNumb; ++i ) { out[i] = rflyColl[i].geti(whichD); }
  return out;
}

/// func
/// ----------------------------------------------------------------------------
/// distribution function : set rdistFunc
/// look at RBeam::rf1
///   - norfac | normalisation factor (fe: to superimpose with rdistHist)
void RBeam::func(double norfac)
{
  /// function name
  string const df[6] = { rname+"_f_x",  rname+"_f_y",  rname+"_f_z", 
                         rname+"_f_px", rname+"_f_py", rname+"_f_pz" }; 
  /// variable name
  string const va[6] = { "x", "y", "z", "p_{x}", "p_{y}", "p_{z}" };
  /// parameters statement 
  double* par = new double[3]; /// functions parameters
  /// main loop
  for( int i=0; i<6; ++i ) { 
    par[0] = norfac; 
    par[1]=rmean[i]; 
    par[2]=rsize[i];
    rfunc[i] = RMath::rf1(df[i], rtype[i], runit[i], par); 
    string xa = va[i] + "  " + "[" + runit[i].getSymbol() + "]";
    rfunc[i].GetXaxis()->SetTitle(xa.c_str());
    rfunc[i].GetXaxis()->CenterTitle();
    /// debug
    //cout << "debug: rfuncxaxis = " << rfunc[i].GetXaxis()->GetTitle() << endl;
  }
  /// release memory
  delete par; par=0;
}

/// funcNH
/// ----------------------------------------------------------------------------
/// normalise rdistFunc to rdistHist.Integral("width")
void RBeam::funcNH()
{
  for( int i=0; i<6; ++i ) {
    RMath::fhscale(rfunc[i], rhist[i]);
  }
}

/// hist
/// ----------------------------------------------------------------------------
/// Sets rhist
/// look RMath::rh1
///   - nbin | bin number in histogram
void RBeam::hist(int const nbin)
{
  /// histogram names
  string hn[6] = { rname+"_h_x",  rname+"_h_y",   rname+"_h_z", 
                   rname+"_h_px", rname+"_h_py",  rname+"_h_pz" }; 
  /// variable name
  string const va[6] = { "x", "y", "z", "p_{x}", "p_{y}", "p_{z}" };
  /// main loop
  for( int i=0; i<6; ++i) {
    double* dcp  = this->flyColl1D(i); /// new double[];
    double  xup  = rfunc[i].GetXmax();
    double  xlow = rfunc[i].GetXmin();
    //for( int i=0; i<rflyNumb; ++i ) { dcp[i] *= runits[i]; }  /// units
    rhist[i] = RMath::rh1( hn[i], dcp, rflyNumb, nbin, runit[i], xlow, xup );
    string xa = va[i] + "  " + "[" + runit[i].getSymbol() + "]";
    rhist[i].GetXaxis()->SetTitle(xa.c_str());
    rhist[i].GetXaxis()->CenterTitle();
    delete[] dcp; dcp=0;
  }
}

/// sizeEmit
/// ----------------------------------------------------------------------------
/// Evaluate 6D beam size from machine parameters
/// For example : sizeMac(x) = sqrt(emittance_x * twiss_beta_x).
///   macParam | pointer to machine parameters:
///     -  par[0] = twiss beta x
///     -  par[1] = twiss beta y
///     -  par[2] = ??? longitudinal parameter that define bunch size
void RBeam::sizeFromEmit(double const mvRef, double const macParam[3])
{
  double* ptr = sizeMac(mvRef, remit, macParam); /// uses new[]
  for(int i=0; i<6; ++i) { 
    rsize[i] = ptr[i]; 
    /// debug 
    /// cout << "db " << i << "  " << "rsize[i] ";
    /// cout << rsize[i]/runit[i] << runit[i].getSymbol() << endl;
  }
  /// release memory
  delete[] ptr; ptr=0;
}







/// static methods  
/// ****************************************************************************

/// static | sizeMac
/// ----------------------------------------------------------------------------
/// computes size from reference momentum, emittance and machine parameters
/// returns a pointer to a 6D size
///   mvRef | reference momentum
///   emitt | emittance (ex, ey, ez)
///   macPa | machine parameters (twissBetax/y/???)
double* RBeam::sizeMac(double const mvRef, double const emitt[3], 
 	                     double const macPa[3])
{
  double* size = new double[6];
  size[0] = sqrt(emitt[0]*macPa[0]);    /// x max
  size[1] = sqrt(emitt[1]*macPa[1]);    /// y max
  //size[2] = sqrt(emitt[2]*macPa[2]);    /// z max ???
  size[2] = 10.*RUnits::mm;    /// force z max
  size[3] = sqrt(emitt[0]/macPa[0]) * mvRef;    /// px max
  size[4] = sqrt(emitt[1]/macPa[1]) * mvRef;    /// py max
  //size[5] = sqrt(emitt[2]/macPa[2]) * mvRef;    /// pz max ???
  size[5] = mvRef;    /// force pz max 
  return size;
}












