# SIRIS4 GEOMETRIC OPTICS WITH DIFFUSE SCATTERERS FRAMEWORK

By:
Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
Department of Physics, University of Helsinki, Finland

## Main programs

- singleparticle Geometrical optics with diffuse scatterers computations for a single Gaussian-random-sphere particle.
- GS Gaussian spheres for creating the discretized meshes for Gaussian sphere particles.

## Compiling, general

The framework is divided into Fortran module files that can be gathered into a siri4-framework subroutine library. The main programs can be compiled and linked using that library. The subroutine library and the most of the main programs are written in Fortran, complying with the 2008 language standard.

There will be a Makefile provided. In general, with GNU gfortran, the following options will are needed/recommended when compiling:

	-ffree-form -std=f2008

### singleparticle Geometric optics computations

SIRIS4 Geometric optics with diffuse scatterers computations for a single Gaussian-random-sphere particle. The SIRIS4-version includes the treatment of inhomogeneous waves in absorbing media. Produces the scattering matrix and efficiences for the particle or the averaged versions over sample of these particles.

If pubishing or distributing results that use this code, please reference both to: Muinonen K, Nousiainen T, Lindqvist H, Muñoz O, and Videen G (2009). Light scattering by Gaussian particles with internal inclusions and roughened surfaces using ray optics. Journal of Quantitative Spectroscopy & Radiative Transfer 110, 1628–1639. DOI:[10.1016/j.jqsrt.2009.03.012](https://doi.org/10.1016/j.jqsrt.2009.03.012), and to Lindqvist H, Martikainen J, Räbinä J, Penttilä A, and Muinonen K (2018). Ray optics in absorbing media with application to ice crystals at near-infrared wavelengths. Journal of Quantitative Spectroscopy & Radiative Transfer 217, 329–337. DOI:[10.1016/j.jqsrt.2018.06.005](https://doi.org/10.1016/j.jqsrt.2018.06.005).

#### Compiling

Option 1: use Makefile

	make singleparticle

Option 2: gfortran compiling with one command:

	gfortran -ffree-form -std=f2008 -o single-particle/siris1p src/sirisconstants.f src/sirismath.f src/sirismaterial.f src/sirisnumint.f src/sirisgeometry.f src/sirisgaussiansphere.f src/sirisradtrans.f src/sirisray.f single-particle/single-particle-main.f

#### Usage

Run from command line. Give the name of the input file as command line argument, e.g.,

	cd single-particle
	./siris1p single-particle-input.in [diffuse-scatterer-scattering-matrix-input-name]

The different parameter options are commented in the example input file.

### GS Gaussian spheres

Produces random Gaussian sphere shapes and discretized mesh representations for them. Can write output in Matlab, IDL IDF, Paraview VTK, and OFF.

If pubishing or distributing results that use this code, please reference to: Muinonen K, Nousiainen T, Fast P, Lumme K, and Peltoniemi JI (1996). Light scattering by Gaussian random particles: Ray optics approximation. Journal of Quantitative Spectroscopy & Radiative Transfer 55(5), 577–601. DOI:[10.1016/0022-4073(96)00003-9](https://doi.org/10.1016/0022-4073(96)00003-9).

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
