### makefile
### 2017/07/10
### Pierre Grandemange
### http://gl.developpez.com/tutoriel/outil/makefile/



### comments
### ----------------------------------------------------------------------------
### 1/	CPLUS_INCLUDE_PATH environnement variable has been set to :
### 			~/programming/repShared/repNamespaces
###				~/programming/repShared/repObjects
### 2/	Carefull on comments, gedit do not see through, see documentation


### particular targets
### ----------------------------------------------------------------------------
.PHONY : clean, mrproper

### variables
### ----------------------------------------------------------------------------
CXX = g++
CXXFLAGS =
LDFLAGS = `root-config --cflags --glibs`
EXEC = xExec.rep
### include folders	: local, namespace, object
LOINC = ../source/include
NSINC = $(repNamespaces)/include
OBINC = $(repObjects)/include
### src folders	: local, namespace, object, data
LOSRC = ../source/src
NSSRC = $(repNamespaces)/src
OBSRC = $(repObjects)/src


## specific!
### --------
### eCooling
ECOOL =$(repData)/coolingData
### environement variable 
export CPLUS_INCLUDE_PATH := ../source/include:$(CPLUS_INCLUDE_PATH):$(ECOOL)


### files list
HEA = $(wildcard $(LOINC)/*.h) \
			$(wildcard $(NSINC)/*.h) $(wildcard $(OBINC)/*.h) \
			$(wildcard $(ECOOL)/*.h)

SRC = ../source/main.cc $(wildcard $(LOSRC)/*.cc) \
      $(wildcard $(NSSRC)/*.cc) $(wildcard $(OBSRC)/*.cc) \
      $(wildcard $(ECOOL)/*.cc)

OBJ = $(SRC:.cc=.o)

### dump some stuff (debug)
### ----------------------------------------------------------------------------
#$(info  $$LOINC is [${LOINC}])
#$(info  $$SRC is [${SRC}])
#$(info  $$HEA is [${HEA}])
#$(info  $$CPLUS_INCLUDE_PATH is [${CPLUS_INCLUDE_PATH}])


### make command: all
### ----------------------------------------------------------------------------
all: $(EXEC) 

### link all .o files
$(EXEC): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

### .o files dependancies
### ---------------------
### main.cc need func.h
### func.h need all object and namespace haders, easier
#main.o: main.cc $(LOINC)/func.h
#func.o: func.cc $(HEA)

#main.o: main.cc
#func.o: func.cc

### .cc -> .o all file auto
%.o: %.cc
	$(CXX) -o $@ -c $< $(LDFLAGS)

	
### make command: clean
### ----------------------------------------------------------------------------
clean:
	rm -rf ../source/main.o
	rm -rf $(LOSRC)/*.o
	rm -rf $(NSSRC)/*.o
	rm -rf $(OBSRC)/*.o
	rm -rf $(ECOOL)/*.o
	
### make command: mrproper
### ----------------------------------------------------------------------------
mrproper: clean
	rm -rf $(EXEC)



