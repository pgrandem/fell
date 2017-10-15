/// RMath.cc
/// ----------------------------------------------------------------------------
/// rep RMath c++ namespace
/// 04/07/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <string>
#include <iostream>
/// root librairies
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
/// rep classes
#include "RUnit.h"
/// rep namespaces
#include "RMath.h"
using namespace std;

namespace RMath {

/// interpolate
/// ----------------------------------------------------------------------------
/// interpolate with function...
double interpol(double first, double last, int n, int i, 
                       string const& model)
{
  double res(-1.);
  if( model=="linear" ) {
    if( n>1 ) {
      double ramp( (last-first)/(n-1) );
      res = first + ramp*i;
    }
    else {
      cout << endl << "rep WARNING : <interpol> func: ";
      cout << "<n> can't be less or equal 1 " << endl; 
    }
  }
  else { 
    cout << endl << "rep WARNING : <interpol> func : model unknowm" << endl; 
  }
  return res;
}

/// y2x()
/// ----------------------------------------------------------------------------
/// compute y from x with model and xi, xf, yi, yf
/// -> linear : y = (yf - yi)/(xf - xi) * (x-xi) + yi
double y2x(double xi, double xf, double x, 
           double yi, double yf, string const& model)
{
  double y(-1);
  if( model=="linear" ) {
    y = (yf - yi)/(xf - xi) * (x-xi) + yi;
  }
  return y;
}



/// gauss - 2017/09/21
/// ----------------------------------------------------------------------------
/// gauss distribution
///   <x>   | the variable
///   <par> | 3 values:
///		  - par[0] : distribution integral value
///		  - par[1] : mean
/// 	  - par[2] : sigma
double gauss(double* x, double* par)
{ return par[0]/(par[2]*sqrt(2*M_PI))*exp(-1./2.*pow((x[0]-par[1])/par[2],2)); }

/// gaussHalf - 2017/09/21
/// ----------------------------------------------------------------------------
/// gauss distribution
///   <x>   | the variable
///   <par> | 3 values:
///		  - par[0] : x0, "mirror edge" of distribution
///		  - par[1] : sigma
/// 	  - par[2] : norm factor (area of gauss = par[2])
double gaussHalf(double* x, double* par)
{ 
  double dist;
  if ( x[0] <= par[1] ) {
    dist = 2. * par[0]/(par[2] * sqrt(2*M_PI)) 
              * exp(-1./2.*pow((x[0]-par[1])/par[2],2));
  } else {dist = 0.; }
  return dist;
}

/// flat - 2017/09/21
/// ----------------------------------------------------------------------------
/// flat distribution
///   <x>   | the variable
///   <par> | two values:
///		  - par[0] : distribution integral value
///		  - par[1] : center of distribution (mean here)
/// 	  - par[2] : size of the flat plateau
double flat(double* x, double* par)
{ 
  double dist(0);
  if( x[0] >= par[1]-par[2]/2. && x[0] <= par[1]+par[2]/2. ) {
  	dist = 1./par[2]*par[0];
  } else { dist = 0.; }
  return dist;
}








/// mean, standartd dev... over many values pointed
/// ****************************************************************************


/// mean
/// ----------------------------------------------------------------------------
/// returns mean value of the x[i]
///   - x | ptr to x[i] values
///   - n | number of values ( 0 <= i < n )
double mean(double* x, int n)
{
  double out(0.);       /// output
  for( int i=0; i<n; ++i ) { out += x[i]; }
  if( n==0 ) { n=1; }
  return (double)out/n;
}

/// sigma
/// ----------------------------------------------------------------------------
/// returns standard deviation value of samples x[i].
///   - x | ptr to x[i] values
///   - n | number of values ( 0 <= i < n )
///   - u | mean value 
double sigma(double* x, int n, double u)
{
  double out(0.); /// output
  for( int i=0; i<n; ++i ) { out += pow( (x[i]-u), 2); }
  if( n==0 ) { n=1; }
  return (double)sqrt(1./n * out);
}

/// mpms
/// ----------------------------------------------------------------------------
/// mean plus/minus sigma
/// returns mean and sigma of a double collection x[i] of size n
///   - x | ptr to x[i] values
///   - n | number of values ( 0 <= i < n )
double* mpms(double* x, int n)
{
  double* out = new double[2];
  out[0] = mean(x, n);
  out[1] = sigma(x, n, out[0]);
  return out;
}






/// rep TF1 from root TF1 functions
/// ****************************************************************************

/// rf1 - 21/09/2017
/// --------------------------------------------------------------------------
/// rep functions 1D. Predefined root TF1 object
/// returns a root TF1 object
///   - name | TF1 name
///   - dtyp | distribution type (gauss, flat, ...)
///   - xuni | x axis unit
///   - para | function parameters. For example, gauss:
///             - par[0] = integral of distribution function
///             - par[1] = center/mean
///             - par[2] = size (WARNING size != sigma)
TF1 rf1(string const& name, string const& dtyp, RUnit const& xunit, double* par)
{
  TF1 out; /// output 
	/// gaussian distribution : size = 4*sigma = xmax
  if( dtyp=="gauss" ) { 
	  double inte = par[0];
	  double mean = par[1]    /xunit;
	  double sigm = par[2]/4. /xunit;
	  double xMax = mean + 4.*sigm;
	  double xmin = mean - 4.*sigm;
	  out = TF1(name.c_str(), RMath::gauss, xmin, xMax, 3);
		out.SetParameter(0, inte);  /// integral value
		out.SetParameter(1, mean);  /// mean/center/y(x)=max...
		out.SetParameter(2, sigm);  /// sigma
	}
	/// flat distribution
	else if( dtyp=="flat" ) { /// flat distribution
	  double inte = par[0];
	  double mean = par[1]  /xunit;
	  double size = par[2]  /xunit;
	  double xMax = mean + size;
	  double xmin = mean - size;
	  out = TF1(name.c_str(), RMath::flat, xmin, xMax, 3);
		out.SetParameter(0, inte);  /// integral value
		out.SetParameter(1, mean);  /// mean/center/y(x)=max...
		out.SetParameter(2, size);  /// size
	}
	/// if "dist" corresponds to nothing, flat distribution 
	else { 
		double inte = 1.; /// integral of dist
    double mean = 0.  /xunit; /// mean/center...
    double size = 1.  /xunit; /// caracteristic "size" of distribution
    double xMax = mean + 1.*size; /// maximum x for the TFunction computation
    double xmin = mean - 1.*size; /// minimum x for the TFunction computation
    out = TF1("noDist", RMath::flat, xmin, xMax, 3);
		out.SetParameter(0, par[0]);  /// integral value
		out.SetParameter(1, par[1]);  /// mean/center/y(x)=max...
		out.SetParameter(2, par[2]);  /// size
	}
	/// TF1 parameters
	out.SetNpx(400);    /// is the output pdf size < 20ko ?!!!
	/// the 3 next lines are commented because axis titles do not remain... 
	//string xaxis = "[" + xunit.getSymbol() + "]";
  //out.GetXaxis()->SetTitle(xaxis.c_str());
  //out.GetXaxis()->CenterTitle();
	/// debug
	//cout << "debug: outxaxis = " << out.GetXaxis()->GetTitle() << endl;
  return out;
}
  
void fhscale(TF1& fun, TH1D const& his)
{
  string xti = fun.GetXaxis()->GetTitle();
  double hsf = his.Integral("width");
  fun.SetParameter(0, hsf);
  fun.GetXaxis()->SetTitle(xti.c_str());
  fun.GetXaxis()->CenterTitle();
}  



/// reprand - 21/09/2017
/// --------------------------------------------------------------------------
/// returns a random distribution of n samples along distfunc
///   - dfunc | distribution function 
///   - xunit | TF1 x axis unit
///   - nsamp | distribution type (gauss, flat, ...)
///   - rseed | random generator seed, 0 means new every time
double* reprand(TF1* const& dfunc, RUnit const& xunit, int nsamp, int ranSeed)
{
  TRandom3 ran(ranSeed);
  double* out = new double[nsamp];
  for( int i=0; i<nsamp; ++i ) { out[i] = dfunc->GetRandom()*xunit; }
  return out;
}

/// rh1 - 21/09/2017
/// --------------------------------------------------------------------------
/// returns a TH1D of nbins filled with a double collection
///   - hn | histo name
///   - dc | double collection values (x[i])
///   - ns | double collection num of samples (0<i=<ns)
///   - nb | number of bins
///   - xu | histo x axis unit
///   - xm | histo x axis min
///   - xM | histo x axis max
TH1D rh1(string const& hn, double* const& dc, int ns, int nb,
         RUnit const& xu, double xm, double xM)
{
  TH1D out(hn.c_str(), hn.c_str(), nb, xm, xM);
  for( int i=0; i<ns; ++i ) { out.Fill( dc[i]/xu ); }
  string xaxis = "[" + xu.getSymbol() + "]";
  out.GetXaxis()->SetTitle(xaxis.c_str());
  return out;
}
















} /// end of namespace


















