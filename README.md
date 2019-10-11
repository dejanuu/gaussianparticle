# SIRIS4 GEOMETRIC OPTICS WITH DIFFUSE SCATTERERS FRAMEWORK

By:
Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
Department of Physics, University of Helsinki, Finland

## Main programs

- GS Gaussian spheres for creating the discretized meshes for Gaussian sphere particles.

## Compiling, general

The framework is divided into Fortran module files that can be gathered into a siri4-framework subroutine library. The main programs can be compiled and linked using that library. The subroutine library and the most of the main programs are written in Fortran, complying with the 2008 language standard.

There will be a Makefile provided. In general, with GNU gfortran, the following options will are needed/recommended when compiling:

	-ffree-form -std=f2008

### GS Gaussian spheres

Produces random Gaussian sphere shapes and discretized mesh representations for them. Can write output in Matlab, IDL IDF, Paraview VTK, and OFF.

If pubishing or distributing results that us this code, please reference to: Muinonen K, Nousiainen T, Fast P, Lumme K, and Peltoniemi JI (1996). Light scattering by Gaussian random particles: Ray optics approximation. Journal of Quantitative Spectroscopy & Radiative Trasnfer 55(5), 577–601. DOI:10.1016/0022-4073(96)00003-9.

#### Compiling

Option 1: use Makefile

	make GS

Option 2: gfortran compiling with one command:

	gfortran -ffree-form -std=f2008 -o GS/GS src/sirisconstants.f src/sirismath.f src/sirisgeometry.f src/sirisgaussiansphere.f GS/GS-main.f

#### Usage

Run from command line. Give the name of the input file as command line argument, e.g.,

	cd GS
	./GS GS-input.in

The different parameter options are commented in the example input file.


