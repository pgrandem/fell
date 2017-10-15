/// prBIDays/readme.txt
/// ----------------------------------------------------------------------------
/// "notebook" like file
/// Pierre Grandemange
/// 09.06.2017
/// ----------------------------------------------------------------------------


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



12.06.2017
/// ----------------------------------------------------------------------------
Reimplementation of my objects
  - RObject : tested. repOject.xx -> RObject.xx



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
  


14.06.2017
/// ----------------------------------------------------------------------------
  -> An easy way to use the correction factor for molecules pressure is to 
  include it in the Molecule object;
  


20.06.2017
/// ----------------------------------------------------------------------------
  Yesterday, start implementing RFlyer and RMachine objects.


