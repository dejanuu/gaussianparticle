# SIRIS4-FRAMEWORK Makefile

# Compiler
COMP ?= gfortran
# Required options for gfortran
FOPT ?= -ffree-form -std=f2008
# Optional choices
# For developing
FOPT += -fimplicit-none -fcheck=bounds,pointer
# Minor optimization
#FOPT += -O1


###############################################################################

LIBFILENAMES = sirisconstants sirismath sirisgeometry sirisgaussiansphere


###############################################################################

all : GS
.PHONY : all clean

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
