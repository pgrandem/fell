### shCompile.sh
### ----------------------------------------------------------------------------
### Script to compile share codes
### Pierre grandemange
### 01/06/2015
### ----------------------------------------------------------------------------



### export compilation needed variables, like folder path
### ----------------------------------------------------------------------------
#export shareProg=~/private/Programming/root/share/prog



### g++ compilation command 
### ----------------------------------------------------------------------------
g++ -I $repObjects/include -I $repNamespaces/include \
    -o xExec.rep main.cc \
    testFunc.cc                               \
    localFunctions.cc                         \
    $repNamespaces/src/RDump.cc               \
    $repNamespaces/src/RParticleList.cc       \
    $repObjects/src/RAtom.cc                  \
    $repObjects/src/RFlyer.cc                 \
    $repObjects/src/RGas.cc                   \
    $repObjects/src/RMachine.cc               \
    $repObjects/src/RMolecule.cc              \
    $repObjects/src/RObject.cc                \
    $repObjects/src/RParticle.cc              \
    `root-config --cflags --glibs`

### g++ simple compilation test command
### ----------------------------------------------------------------------------
#g++ -o xExec.rep main.cc `root-config --cflags --glibs`



### trash
### ----------------------------------------------------------------------------
### $repObjects/src/RInteraction_FlyerGas.cc  \
    
