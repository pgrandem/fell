/// RObject.h
/// ----------------------------------------------------------------------------
/// rep RObject c++ object
/// 12/06/2017
/// Pierre Grandemange
/// ----------------------------------------------------------------------------

#ifndef DEF_ROBJECT
#define DEF_ROBJECT


/// includes
/// ----------------------------------------------------------------------------
/// standard libraries
#include <iostream>
#include <string>
/// root librairies
/// rep functions
/// rep objects
/// local functions



class RObject
{
  public :
  /// constructor, destructor, copy 
  /// **************************************************************************
  RObject();                          /// constructor
  RObject(const std::string& name);   /// constructor
  virtual ~RObject();                 /// destructor
  RObject(RObject const& objToCopy);  /// copy constructor
    
  /// accessors
  /// **************************************************************************
  void          setName(const std::string& name) { rname = name; };
  std::string   getName() const                  { return rname; };
    
  /// methods
  /// **************************************************************************
  /// dump properties
  /// --------------------------------------------------------------------------
  virtual void dump(std::ostream &flux=std::cout) const;
  
  
  
  /// operators
  /// **************************************************************************
  /// sort all Rep object by name!
  /// implemented to sort maps of machine elements
  RObject&  operator=(const RObject& otherObj);
  bool      operator<(const RObject& otherObj) const; /// use string < operator
  
  /// static methods  
  /// **************************************************************************
  static int  nb();     /// return the number of instance created
  
  /// attribute
  /// **************************************************************************
  protected :
  std::string rname;     /// [ ] object name
  
  /// static attributes  
  /// **************************************************************************
  static int  rnb;   /// number of running instances

};



/// external operators
/// ****************************************************************************
std::ostream& operator<<(std::ostream &flux, RObject const& object);

#endif /// DEF_ROBJECT



