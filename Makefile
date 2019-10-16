# SIRIS4-FRAMEWORK Makefile

# Compiler
COMP ?= gfortran
# Required options for gfortran
FOPT = -ffree-form -std=f2008
# Optional choices
# For developing
#FOPT += -fimplicit-none -fcheck=bounds,pointer -Wall -Wno-maybe-uninitialized
# Minor optimization
#FOPT += -O1
# Major optimization
FOPT += -Ofast -mtune=native
# For debugging
#FOPT += -g -ffpe-trap=invalid,zero,overflow,underflow

###############################################################################

MKDIR = mkdir -p
DIRECTORIES = lib mod
LIBFILENAMES = sirisconstants sirismath sirismaterial sirisnumint sirisgeometry sirisgaussiansphere sirisray sirisradtrans

###############################################################################

all : directories GS singleparticle
.PHONY : all clean directories

singleparticle : single-particle/single-particle-main.f siris4lib
	$(COMP) $(FOPT) -o single-particle/singleparticle -J mod -L lib single-particle/single-particle-main.f -l siris4 

GS : GS/GS-main.f siris4lib
	$(COMP) $(FOPT) -o GS/GS -J mod -L lib GS/GS-main.f -l siris4 

siris4lib : lib/libsiris4.a
	
lib/libsiris4.a : $(addprefix lib/,$(addsuffix .o,$(LIBFILENAMES)))
	ar crv -s lib/libsiris4.a $^

lib/%.o : src/%.f
	$(COMP) $(FOPT) -c -J mod -o lib/$*.o src/$*.f

clean : 
	rm -rf lib/*.o
	rm -rf lib/*.a
	rm -rf mod/*.mod
	rm -rf GS/GS.exe
	rm -rf GS/GS
	rm -rf single-particle/singleparticle.exe
	rm -rf single-particle/singleparticle

directories: ${DIRECTORIES}

${DIRECTORIES}:
	${MKDIR} ${DIRECTORIES}
