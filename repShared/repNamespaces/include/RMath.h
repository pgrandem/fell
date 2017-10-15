/// RMath.h
/// ----------------------------------------------------------------------------
/// rep RMath c++ namespace
/// 03/07/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------


#ifndef DEF_RMATH
#define DEF_RMATH


/// includes
/// ----------------------------------------------------------------------------
/// standard libraries
#include <cmath>
/// root librairies
#include "TF1.h"
#include "TH1D.h"
/// rep functions
/// rep objects
#include "RUnit.h"
/// local functions

namespace RMath {

/// interpol()
/// ----------------------------------------------------------------------------
/// interpolate with model
double interpol(double first, double last, int npts, int iterator, 
                std::string const& model="linear");
/// y2x()
/// ----------------------------------------------------------------------------
/// compute y from x with model and xi, xf, yi, yf
/// -> linear : y = (yf - yi)/(xf - xi) * (x-xi) + yi
double y2x(double xi, double xf, double x, 
           double yi, double yf, std::string const& model="linear");




/// c type functions distributions (syntax found for <ROOT TF1>)
/// ****************************************************************************

/// gauss - 2017/09/21
/// ----------------------------------------------------------------------------
/// gauss distribution
///   <x>   | the variable
///   <par> | 3 values:
///		  - par[0] : distribution integral value
///		  - par[1] : center of distribution (mean here)
/// 	  - par[2] : sigma
double gauss(double* x, double* par);

/// gaussHalf - 2017/09/21
/// ----------------------------------------------------------------------------
/// gauss distribution
///   <x>   | the variable
///   <par> | 3 values:
///		  - par[0] : distribution integral value
///		  - par[1] : edge(center of normal dist) value
/// 	  - par[2] : sigma
double gaussHalf(double* x, double* par);

/// flat - 2017/09/21
/// ----------------------------------------------------------------------------
/// flat distribution
///   <x>   | the variable
///   <par> | 3 values:
///		  - par[0] : distribution integral value
///		  - par[1] : center of distribution (mean here)
/// 	  - par[2] : sigma
double flat(double* x, double* par);





/// mean, standartd dev... over many values pointed
/// ****************************************************************************

/// mean
/// ----------------------------------------------------------------------------
/// returns mean value of the x[i]
///   - x | ptr to x[i] values
///   - n | number of values ( 0 <= i < n )
double mean(double* x, int n);

/// sigma
/// ----------------------------------------------------------------------------
/// returns standard deviation value of samples x[i].
///   - x | ptr to x[i] values
///   - n | number of values ( 0 <= i < n )
///   - u | mean value 
double sigma(double* x, int n, double u);

/// mpms
/// ----------------------------------------------------------------------------
/// mean plus/minus sigma
/// returns mean and sigma of a double collection x[i] of size n
///   - x | ptr to x[i] values
///   - n | number of values ( 0 <= i < n )
double* mpms(double* x, int n);



/// rep "root predefined" object
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
///             - par[2] = size of distribution (WARNING size != sigma)
TF1 rf1(std::string const& name, std::string const& dtyp, 
        RUnit const& xunit, double* par);

void fhscale(TF1& fun, TH1D const& his);


/// reprand - 21/09/2017
/// --------------------------------------------------------------------------
/// returns a random distribution of n samples along distfunc
///   - dfunc | distribution function 
///   - xunit | TF1 x axis unit
///   - nsamp | distribution type (gauss, flat, ...)
///   - rseed | random generator seed, 0 means new every time
double* reprand(TF1* const& dfunc, RUnit const& xunit, int nsamp, 
                int ranSeed=0);

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
TH1D rh1(std::string const& hn, double* const& dc, int ns, int nb, 
         RUnit const& xu, double xm, double xM);



} /// end of namespace 

#endif /// RMATH
