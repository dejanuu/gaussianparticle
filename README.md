# SIRIS4 GEOMETRIC OPTICS WITH DIFFUSE SCATTERERS FRAMEWORK

By:
Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
Department of Physics, University of Helsinki, Finland

## Main programs

- GS Gaussian spheres for creating the discretized meshes for Gaussian sphere particles.

## Compiling, general

The framework is divided into Foetran module files that can be gathered into a siri4-framework subroutine library. The main programs can be compiled and linked using that library. The subroutine library and the most of the main programs are written in Fortran, complying with 2008 language standard.

There will be a Makefile provided. In general, with GNU gfortran, the following options will are needed/recommended when compiling: '-ffree-form -std=f2008 -fimplicit-none'.

### GS Gaussian spheres

Produces random Gaussian sphere shapes and discretized mesh representations for them. Can write output in Matlab, IDL IDF, Paraview VTK, and OFF.

#### Compiling

Option 1: use Makefile.

Option 2: gfortran compiling with one command:
	 gfortran -ffree-form -std=f2008 -fimplicit-none -o GS sirisconstants.f sirismath.f sirisgeometry.f sirisgaussiansphere.f GS-main.f

#### Usage

XXX


