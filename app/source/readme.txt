/// prBIDays/readme.txt
/// ----------------------------------------------------------------------------
/// "notebook" like file
/// Pierre Grandemange
/// 09.06.2017
/// ----------------------------------------------------------------------------




21/09/2017
/// ----------------------------------------------------------------------------


20/09/2017
/// ----------------------------------------------------------------------------
"Working/compiling" new version of beam avalaible.
Changing RGas class to adapt.


19/09/2017
/// ----------------------------------------------------------------------------
readme.txt (this file) moved from build to source folder.
RBeam class changed a lot.
  -> First use in <testBeamDist()> test function.
  -> Still some old method used in RGas class that must be fixed.
Integration of millisecond counter in main and sub function routines.


17/09/2017
/// ----------------------------------------------------------------------------
New backup :
  -> /home/rep/Documents/_backup/fell_source_20170917
  -> /home/rep/Documents/_backup/repObjects_20170917
  


16/09/2017
/// ----------------------------------------------------------------------------
The soft is now working on my personal machine (rep-UX390UAK).
The development was stopped during "beam" class implementation:
 -> what's a beam? a profile? a distribution? ...
 
16h14.This is where is my reflection:
  Flyer: A flyer is a cinetic particle defined by [x,y,z,px,py,pz].
    -> If defined by rel. param. beta, then [0, 0, 0, 0, 0, gamma*beta*m0*c]. 
  Beam: defined by emittance/array of flyers and particles. 
    -> Allows to create distributions of flyers.

20h51
RFlyer has been modified deeply.
RBeam is next, I started messing with it. It means problems with RGas class.
  -> To be checked and reimplemented.
<make mrproper> realised and <make> again. It's compiling. 
May(99% sure) introduced bugs in former <test>Functions.



12/09/2017
/// ----------------------------------------------------------------------------
Migration of this software to my personnal computer rep-UX390UAK
today's goal : make this software to work on this machine.



20.06.2017
/// ----------------------------------------------------------------------------
Yesterday, start implementing RFlyer and RMachine objects.



15.06.2017
/// ----------------------------------------------------------------------------
  -> An easy way to use the correction factor for molecules pressure is to 
  include it in the Molecule object;
  


14.06.2017
/// ----------------------------------------------------------------------------
Reimplementation of objects
  - RObject
  - RParticle
  - RAtom
Implementation of namespaces
  - RUnits.h
  - RParticleList.h
  - RParticleList.cc

Explanations about relative atomic mass/standard atomic weight "assimilated" ?!
  ->  http://physics.nist.gov/cgi-bin/Compositions/
      stand_alone.pl?ele=&ascii=html&isotype=some

Reminder : old objects to be modified can be find in:
  repShared/_save/repObjects_20170609

Molecular mass values:
  http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440633&Units=SI&Mask=21
  


12.06.2017
/// ----------------------------------------------------------------------------
Reimplementation of my objects
  - RObject : tested. repOject.xx -> RObject.xx



09.06.2017
/// ----------------------------------------------------------------------------
En essayant de reprendre le programme BIPM pour BIDays:
  - pb de compilation
  - j'essaie de faire fusionner mes objets entre:
    - prBIPM
    - prBeams
    -geant4/TestEm5/detRes -> le plus aboutit
    
  - les problemes de compilation incluent notamment des boucles de declaration
    d'objets. J'apprends à utiliser les "forward declaration"...

Les problèmes de déclaration sont tellement lourds, et les cross references
  tellement nombreuses qu'il vaut mieux reprendre doucement...

-> création de ce dossier!



