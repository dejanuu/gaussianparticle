# SIRIS4-FRAMEWORK Makefile

# Compiler
COMP ?= gfortran
CPPCOMP ?= g++
# Required options for gfortran
FOPT = -ffree-form -std=f2008
# Optional choices
# For developing
#FOPT += -fimplicit-none -fcheck=bounds,pointer -Wall -Wno-maybe-uninitialized
# Minor optimization
FOPT += -O1
# Major optimization
#FOPT += -Ofast -mtune=native -ffpe-summary=none
# For debugging
#FOPT += -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow

###############################################################################

MKDIR = mkdir -p
DIRECTORIES = lib mod
LIBFILENAMES = sirisconstants sirismath sirismaterial sirisnumint sirisgeometry sirisgaussiansphere sirisray sirisradtrans
MPFFILES = mathroutines splinetools sirisinterface
MPCFILES = detector geometry inputreader materials meshreader mray offreader outputwriter physicsengine rng smaterial sray tracer

###############################################################################

all : directories GS singleparticle singletwolayer
.PHONY : all clean directories

singletwolayer : directories single-two-layer-particle/single-two-layer-main.f siris4lib
	$(COMP) $(FOPT) -o single-two-layer-particle/siris2l -J mod -L lib single-two-layer-particle/single-two-layer-main.f -lsiris4 

singleparticle : directories single-particle/single-particle-main.f siris4lib
	$(COMP) $(FOPT) -o single-particle/siris1p -J mod -L lib single-particle/single-particle-main.f -lsiris4 

GS : directories GS/GS-main.f siris4lib
	$(COMP) $(FOPT) -o GS/GS -J mod -L lib GS/GS-main.f -lsiris4

multiparticle : directories multi-particle/main.cpp $(addprefix multi-particle/,$(addsuffix .o,$(MPFFILES))) $(addprefix multi-particle/,$(addsuffix .o,$(MPCFILES))) siris4lib
	$(CPPCOMP) -o multi-particle/sirismp -L lib $(addprefix multi-particle/,$(addsuffix .o,$(MPFFILES))) $(addprefix multi-particle/,$(addsuffix .o,$(MPCFILES))) multi-particle/main.cpp -lCGAL -lsiris4 -lgfortran 

siris4lib : directories lib/libsiris4.a
	
lib/libsiris4.a : $(addprefix lib/,$(addsuffix .o,$(LIBFILENAMES))) 
	ar crv -s lib/libsiris4.a $^

lib/%.o : src/%.f
	$(COMP) $(FOPT) -c -J mod -o lib/$*.o src/$*.f

multi-particle/%.o : multi-particle/%.f
	$(COMP) $(FOPT) -c -J mod -o multi-particle/$*.o multi-particle/$*.f

multi-particle/%.o : multi-particle/%.cpp
	$(CPPCOMP) -c -o multi-particle/$*.o multi-particle/$*.cpp

clean : 
	rm -rf lib/*.o
	rm -rf lib/*.a
	rm -rf mod/*.mod
	rm -rf GS/GS.exe
	rm -rf GS/GS
	rm -rf single-particle/siris1p.exe
	rm -rf single-particle/siris1p
	rm -rf single-two-layer-particle/siris2l.exe
	rm -rf single-two-layer-particle/siris2l
	rm -rf multi-particle/*.o
	rm -rf multi-particle/sirismp.exe
	rm -rf multi-particle/sirismp

directories: ${DIRECTORIES}

${DIRECTORIES}:
	${MKDIR} ${DIRECTORIES}
