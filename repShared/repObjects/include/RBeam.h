/// RBeam.h
/// ----------------------------------------------------------------------------
/// rep RBeam c++ class 
/// Pierre GRANDEMANGE
/// 13/07/2017
/// ----------------------------------------------------------------------------


#ifndef DEF_RBEAM
#define DEF_RBEAM


/// includes
/// ----------------------------------------------------------------------------
/// standard library
#include <string>
/// root classes
#include "TF1.h"
#include "TH1D.h"
/// rep classes
#include "RObject.h"
#include "RParticle.h"
#include "RUnit.h"

class RFlyer;
class RMachine;
class RGas;
//class RUnit;

class RBeam : public RObject
{
	/// attributes
	/// **************************************************************************
	protected:
	double			remit[3];	  /// [um, um, ???] beam dist physical emittance 3D
	RFlyer*			rflyColl;   /// [ ] pointer to flyer collection (array)
	int    	    rflyNumb;	  /// [ ] beam dist number of flyer
	TF1         rfunc[6];   /// [ ] beam dist function(type, size) 6D
	TH1D        rhist[6];   /// [ ] beam dist histogram(function, nb)
	double      rmean[6];   /// [ ] beam dist 6D  "mean"
 	RParticle*  rparType;   /// [ ] beam dist particle type
	double      rsize[6];   /// [RUnits, RU...] beam dist "size" 6D
	std::string rtype[6];		/// [ ] beam dist type (gaussian, flat, etc...) 6D
	RUnit       runit[6];   /// [ ] units for plotting (func and hist)
	
 	
 	/// constructor, destructor, copy 
	/// **************************************************************************
	public:
	RBeam(std::string const& name, int npb, RParticle* par, 
	      std::string const distType[6],  double mv0, double emit[3], 
	      double macPar[3]);
  
  RBeam(std::string const& name, int npb, RParticle* par, 
	      std::string const distType[6], double mean[6], double size[6]);
  ~RBeam();
	
  
  
  
  /// accessors
  /// **************************************************************************
 	public:
 	TF1* getFunc() const          { return (TF1*)&rfunc[0]; }
 	void setFunc(TF1 const df[6]) { for(int i=0;i<6;++i) rfunc[i]=df[i]; } 
 	
 	TH1D* getHist() const { return (TH1D*)&rhist[0]; }
 	void  setHist( TH1D const h[6] ) 
 	                 { for(int i=0;i<6;++i) rhist[i]=h[i]; } 
 	
 	std::string* getType() const { return (std::string*)&rtype[0]; }
 	void         setType(std::string const dt[6]) 
 	                        { for(int i=0;i<6;++i) rtype[i]=dt[i]; } 
 	
 	double* getEmit() const { return (double*)&remit[0]; }
 	void    setEmit(double emit[3])  
 	                    { for(int i=0;i<3;++i) remit[i] = emit[i]; }
 	
 	RFlyer* getFlyColl()  const             { return rflyColl; }
 	void    setFlyColl( RFlyer* flyerColl ) { rflyColl = flyerColl;  }
 	
 	int   getN() const  { return rflyNumb; }
	void  setN(int n)		{ rflyNumb=n; }
	
	double* getMean() const         { return (double*)&rmean[0]; }
  void    setMean(double mean[6]) { for(int i=0; i<6; ++i) rmean[i]=mean[i]; } 
	
	RParticle*  getParticle() const           { return rparType; }
	void        setParticle( RParticle* par ) { rparType = par; }
	
	double* getSize() const         { return (double*)&rsize[0]; }
  void    setSize(double size[6]) { for(int i=0; i<6; ++i) rsize[i]=size[i]; } 
	
	/// for dvlpt issues to fix with RGas class using former RBeam methods
	RFlyer* getFlyer() { return rflyColl; }
	//double  getNpb()   { return this->getN(); }
	
	
	
	
	/// indirect accessors
	/// --------------------------------------------------------------------------
  const TH1D* hx()   const { return (const TH1D*)&rhist[0]; }
	const TH1D* hy()   const { return (const TH1D*)&rhist[1]; } 
	const TH1D* hz()   const { return (const TH1D*)&rhist[2]; } 
	const TH1D* hpx()  const { return (const TH1D*)&rhist[3]; } 
	const TH1D* hpy()  const { return (const TH1D*)&rhist[4]; } 
	const TH1D* hpz()  const { return (const TH1D*)&rhist[5]; } 
	
 	
 	
 	
 	/// dump methods
  /// **************************************************************************
 	public:
 	/// dump line
  /// --------------------------------------------------------------------------
  /// dump main beam properties on one line
  void dumpLine(std::ostream& flux=std::cout) const;
  
  /// dumpFlyColl
  /// --------------------------------------------------------------------------
  /// dump flyer collection x y z px py pz 
  void dumpFlyColl(int modulo=1, std::ostream& flux=std::cout) const;
  
  
  
 	/// other methods
  /// **************************************************************************
 	public:
	/// flyColl
  /// --------------------------------------------------------------------------
  /// Set rflyColl from rfunc and rflyNumb.
  /// Look at RMath::random.
  void flyColl(int ranSeed=0);
	
	/// flyColl1D
	/// --------------------------------------------------------------------------
  /// returns a double ptr to 1D(over 6) of a flyer collection
	///   -whichD | which dim to return (0 to 5)
	double* flyColl1D( int whichD) const;
  double* x()   const { flyColl1D(0); }
  double* y()   const { flyColl1D(1); }
  double* z()   const { flyColl1D(2); }
  double* px()  const { flyColl1D(3); }
  double* py()  const { flyColl1D(4); }
  double* pz()  const { flyColl1D(5); }
  
	/// func
  /// --------------------------------------------------------------------------
  /// Sets rfunc from rsize, and rtype (and norfac).
  /// Look at RBeam::rf1.
  ///   - norfac | normalisation factor (fe: to superimpose with rdistHist)
  void func(double norfac=1.);
	
	/// funcNH
	/// --------------------------------------------------------------------------
  /// normalise rdistFunc to rdistHist.Integral("width")
	void funcNH(); /// func : normalize o histogram
	
	/// hist
  /// --------------------------------------------------------------------------
  /// Sets rhist
  /// look RMath::rh1
  ///   - nbin | bin number in histogram
  void hist(int const nbin);
	
  /// sizeEmit
  /// --------------------------------------------------------------------------
  /// set rsize from remit and machine parameters
  /// For example : sizeMac(x) = sqrt(emittance_x * twiss_beta_x).
  ///   mvRef | beam reference particle momentum
  ///   macParam | pointer to machine parameters:
  ///     -  par[0] = twiss beta x
  ///     -  par[1] = twiss beta y
  ///     -  par[2] = ??? longitudinal parameter that define bunch size
  void sizeFromEmit(double const mvRef, double const macPa[3]);





	/// static methods  
  /// **************************************************************************
 	public:
 	/// static | sizeMac
 	/// --------------------------------------------------------------------------
  /// computes size from reference momentum, emittance and machine parameters
  /// returns a pointer to a 6D size
  ///   mvRef | reference momentum
  ///   emitt | emittance (ex, ey, ez)
  ///   macPa | machine parameters (twissBetax/y/???)
  static double* sizeMac(double const mvRef, double const emitt[3], 
 	                        double const macPa[3]);




};

#endif /// RBEAM













