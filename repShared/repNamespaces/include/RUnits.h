/// RUnits.h
/// ----------------------------------------------------------------------------
/// rep RUnits c++ namespace
/// 13/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------



/// basic units are SI
/// ----------------------------------------------------------------------------
/// category    base_unit     symbol
/// --------    ---------     ------
/// weight      kilogram      kg
/// length      meter         m
/// time        second        s
/// energy      joule         J
/// charge      coulomb       C
/// ...


#ifndef DEF_RUNITS
#define DEF_RUNITS


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <cmath>
#include "RUnit.h"

namespace RUnits {

/// SI units
/// ----------------------------------------------------------------------------
/// -> https://en.wikipedia.org/wiki/International_System_of_Units
static const double kilogram  = 1.;           /// weight
static const double kg        = 1.*kilogram;
//static const double second    = 1.;           
//static const double s         = 1.*second;
static const RUnit s("second", "s", "time", 1.);  /// time
static const RUnit second = s;
static const RUnit m("meter", "m", "length", 1.); /// length
static const RUnit meter = m;
static const RUnit A("ampere", "A", "current", 1.); /// length
static const RUnit ampere = A;
//static const double meter     = 1.;           
//static const double m         = 1.*meter;
//static const double ampere    = 1.;           /// electrical current
//static const double A         = 1.*ampere;
static const double kelvin    = 1.;           /// thermodynamic temperature
static const double K         = 1.*kelvin;
static const double mole      = 1.;           /// amount of substance
static const double mol       = 1.*mole;
static const double candela   = 1.;           /// luminous intensity
static const double cd        = 1.*candela;



/// fundamental/elementary constants
/// ----------------------------------------------------------------------------
static const double c     = 299792458.*m/s;                 /// speed unit
static const double kB    = 1.38064852e-23*kg*m*m/(s*s)/K;  /// <=> J/K
static const double e     = 1.60217662e-19*A*s;             /// A*s = C
static const double eps0  = 8.85418782e-12 * pow(A.v(),2) * pow(s.v(),4) / kg 
                                           / pow(m.v(),3);
static const double one   = 1.; /// no unit

/// convertors
/// ----------------------------------------------------------------------------
static const double kilo  	= 1.e3;
static const double Mega  	= 1.e6;
static const double percent = 1.e-2;
static const double milli 	= 1.e-3;
static const double micro 	= 1.e-6;
static const double nano 		= 1.e-9;








/// action
/// ----------------------------------------------------------------------------
//static const double

/// charge
/// ----------------------------------------------------------------------------
static const double coulomb = 1.*A*s;
static const double C       = 1.*coulomb;

/// electric potential
/// ----------------------------------------------------------------------------
static const double volt  = 1.*kg*m*m/A/(s*s*s);
static const double V     = 1.*volt;

/// energy
/// ----------------------------------------------------------------------------
static const double joule = 1.*kg*m*m/(s*s);
static const double J     = 1.*joule;
static const double eV    = 1.*e*V;
static const double keV   = 1.*kilo*eV;
static const double MeV   = 1.*Mega*eV;

/// momentum
/// ----------------------------------------------------------------------------
//static const double MeV_c = 1.*MeV/c;
static const RUnit MeV_c("MeV/c", "MeV/c", "momentum", 1.*MeV/c);
static const double keV_c = 1.*keV/c;

/// pressure
/// ----------------------------------------------------------------------------
static const double pascal  = 1.*kg/m/(s*s); /// kg/m/s2
static const double Pa      = 1.*pascal;
static const double bar     = 1.e5*pascal;
static const double mbar    = 1.*milli*bar;

/// length
/// ----------------------------------------------------------------------------
static const double angstrom  	= 1.00001501e-10*meter;
static const double centimeter  = 0.01*meter;
static const double cm  				= 1.*centimeter;
//static const double millimeter  = 1.e-3*meter;
//static const double mm  				= 1.*millimeter;
static const RUnit mm("millimeter", "mm", "length", 1.*milli*m);
//static const double micrometer  = 1.*micro*meter;
//static const double um  				= 1.*micrometer;
static const RUnit um("micrometer", "um", "length", 1.*micro*m);


/// surface
/// ----------------------------------------------------------------------------
static const double square_meter  = 1.*meter*1.*meter;
static const double m2            = 1.*square_meter;

/// temperature
/// ----------------------------------------------------------------------------

/// time
/// ----------------------------------------------------------------------------
static const double millisecond = 1.*milli*second;
static const double ms					= 1.*millisecond;
static const double microsecond = 1.*micro*second;
static const double us          = 1.*microsecond;
static const double nanosecond  = 1.*nano*second;
static const double	ns					= 1.*nanosecond;

/// speed
/// ----------------------------------------------------------------------------
static const double m_s = 1.*m/s;
static const double m_us = 1.*m/us;
static const double um_s = 1.*micro*m_s;

/// volume
/// ----------------------------------------------------------------------------
static const double cube_meter  = 1.*meter*1.*meter*1.*meter;
static const double m3          = 1.*cube_meter;

/// weight
/// ----------------------------------------------------------------------------
static const double MeV_c2  = 1.*MeV/(c*c);
static const double amu     = 1.66053892e-27*kg;

}
#endif /// DEF_RUNITS
