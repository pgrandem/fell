/// RFlyer.h
/// ----------------------------------------------------------------------------
/// rep RFlyer c++ class
/// Pierre Grandemange
/// 19/06/2017
/// ----------------------------------------------------------------------------


#ifndef DEF_RFLYER
#define DEF_RFLYER


/// includes
/// ----------------------------------------------------------------------------
/// standard library
/// root classes
/// rep namespaces
/// rep classes
#include "RParticle.h"
#include "RUnits.h"


class RFlyer : public RParticle
{
  /// attribute
  /// **************************************************************************
  protected:
  double    rpos6D[6];  /// [ ] x, y, z, px, py, pzpublic :
  
  
  
  
  /// constructor, destructor, copy 
  /// **************************************************************************
  public:
  RFlyer();
  RFlyer(std::string const& name, std::string const& symbol, 
         double charge, double mass, double beta);
  RFlyer(RParticle* const& particle, double beta);
  RFlyer(RParticle* const& particle, double const pos6D[6]);
  virtual ~RFlyer();
  
  
  /// accessors
  /// **************************************************************************
  public:
  //double getBeta() const        { return rbeta; };
  //double setBeta(double beta)   { rbeta = beta; };
  double* getPos6D() const { return (double*)&rpos6D[0]; };
  void    setPos6D( double x, double y, double z, 
                    double px, double py, double pz ) 
                    { rpos6D[0]=x;  rpos6D[1]=y;  rpos6D[2]=z; 
                      rpos6D[3]=px; rpos6D[4]=py; rpos6D[5]=pz; };
  
  
  /// dump method
  /// **************************************************************************
  public:
  /// dump
  /// --------------------------------------------------------------------------
  /// dump some properties
  virtual void  dump(     std::ostream &flux=std::cout  )  const;
  void          dumpBeta( std::ostream &flux=std::cout  )  const;
  void          dump6D(   std::ostream &flux=std::cout  )  const;
  void          line6D(   std::ostream &flux=std::cout  )  const;
  
  
  
  
  
  
  
  
  
  /// normal methods
  /// **************************************************************************
  
  public:
  ///particle
  /// --------------------------------------------------------------------------
  RParticle* particle() const;        /// build particle from flyer
  void       particle(RParticle* par);/// set the particle properties of a flyer
  
  double geti(int i)  const { return rpos6D[i]; }
  double x()  const { return rpos6D[0]; }
  double y()  const { return rpos6D[1]; }
  double z()  const { return rpos6D[2]; }
  double px() const { return rpos6D[3]; }
  double py() const { return rpos6D[4]; }
  double pz() const { return rpos6D[5]; }
  double p()  const { return sqrt( pow(px(),2)+pow(py(),2)+pow(pz(),2) ); }
  double xp() const { return px()/p(); }
  double yp() const { return py()/p(); }
  
  double distance() const { return sqrt(pow(x(),2)+ pow(y(),2)+ pow(z(),2)); };
  double momentum() const { return sqrt(pow(px(),2)+pow(py(),2)+pow(pz(),2)); };
  double beta()   const { return momentumToBeta(momentum(), rmass); };
  double gamma()  const { return betaToGamma(beta()); };
  
  double ekin()     const { return betaToEkin(beta(), rmass); };
  double e0()       const { return massToE0(rmass);           };
  double velocity() const { return betaToVelocity(beta());    };
  
  void momentum(double p) { setPos6D( 0., 0., 0., 0., 0., p); };
  void ekin(double ekin)  { setPos6D( 0., 0., 0., 0., 0., 
                            ekinToMomentum(ekin, rmass)   );  };
  
  
  
  
  
  
  
  
  /// static methods  
  /// **************************************************************************
  static double betaToEkin(double beta, double mass);
  static double betaToGamma(double beta);
  static double betaToMomentum(double beta, double mass); 
  static double betaToVelocity(double beta); 
  static double ekinToBeta(double ekin, double mass);
  static double ekinToGamma(double ekin, double mass);
  static double ekinToMomentum(double ekin, double mass);
  static double gammaToBeta(double gamma);
  static double gammaToEkin(double gamma, double mass);
  static double massToE0(double mass);
  static double momentumToBeta(double momentum, double mass);
  static double momentumToEkin(double momentum, double mass);
  
  
 
 
 
};

#endif /// RFlyer









