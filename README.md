# SIRIS4 GEOMETRIC OPTICS WITH DIFFUSE SCATTERERS FRAMEWORK

By:
Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä

Department of Physics, University of Helsinki, Finland

## Main programs

- *siris1p* Geometric optics with diffuse scatterers computations for a single Gaussian-random-sphere particle.
- *siris2l* Geometric optics with diffuse scatterers computations for a single Gaussian-random-sphere particle having core-mantle geometry.
- *sirismp* Geometric optics with diffuse scatterers computations for multiple particles.
- *GS* Gaussian spheres for creating the discretized meshes for Gaussian sphere particles.

## Compiling, general

The framework is divided into Fortran module files that can be gathered into a siri4-framework subroutine library. The main programs can be compiled and linked using that library. The subroutine library and the most of the main programs are written in Fortran, complying with the 2008 language standard.

There will be a Makefile provided. In general, with GNU gfortran, the following options will are needed/recommended when compiling:

	-ffree-form -std=f2008

### siris1p Geometric optics computations

SIRIS4 Geometric optics with diffuse scatterers computations for a single Gaussian-random-sphere particle. The SIRIS4-version includes the treatment of inhomogeneous waves in absorbing media. Produces the scattering matrix and efficiences for the particle or the averaged versions over sample of these particles.

If pubishing or distributing results that use this code, please reference both to \[2\] and \[3\].

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

### siris2l Geometric optics computations for two-layer particle

SIRIS4 Geometric optics with diffuse scatterers computations for a single Gaussian-random-sphere particle with core-mantle structure. The SIRIS4-version includes the treatment of inhomogeneous waves in absorbing media. Produces the scattering matrix and efficiences for the particle or the averaged versions over sample of these particles.

If pubishing or distributing results that use this code, please reference both to \[2\] and \[3\].

#### Compiling

Option 1: use Makefile

	make singletwolayer

Option 2: gfortran compiling with one command:

	gfortran -ffree-form -std=f2008 -o single-two-layer-particle/siris2l src/sirisconstants.f src/sirismath.f src/sirismaterial.f src/sirisnumint.f src/sirisgeometry.f src/sirisgaussiansphere.f src/sirisradtrans.f src/sirisray.f single-two-layer-particle/single-two-layer-main.f

#### Usage

Run from command line. Give the name of the input file as command line argument, e.g.,

	cd single-two-layer-particle
	./siris2l single-two-layer-input.in [mantle-diffuse-scatterer-scattering-matrix-input-name] [core-diffuse-scatterer-scattering-matrix-input-name]

The different parameter options are commented in the example input file.

### sirismp Geometric optics computations for multiple particles

SIRIS4 Geometric optics with diffuse scatterers computations for multiple particles. The SIRIS4-version includes the treatment of inhomogeneous waves in absorbing media. Produces the scattering matrix and efficiences for the collection of particles.

If pubishing or distributing results that use this code, please reference to \[4\].

#### Compiling

Option 1: use Makefile

	make multiparticle

Option 2: use CMake

	cd multi-particle
	mkdir build
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=Release
	make
	cp SIRISn ../sirismp

With both options, you will need the C++-libraries CGAL (https://www.cgal.org/) and Boost (https://www.boost.org/) installed in your system. Note, if compilation problems with include files or the definition of M_PI, edit Makefile and/or multi-particle/definitions.hpp.

#### Usage

Run from command line. Give the name of the input file as command line argument, e.g.,

	cd multi-particle/example
	../sirismp input.in

The possible arguments are listed below, and are given in the file in format:
```
ARGUMENT VALUE
```

Possible arguments are:

```
nrays 1000000
```
Number of rays

```
max_scattering 200
```
Maximum number of scattering events

```
killswitch_start 70
```
Prevent rays splitting to refracted and reflected (only one of these happens) rays after N scattering events 

```
nbins 80
```
Number of theta angle bings

```
nbins_fine_details_start 180
```
This is used to print finer details at the backscattering direction

```
prevent_TR 0
```
Prevent total reflection creating refracted rays

```
check_time_after_nrays 1000
```
Check time after N rays

```
allocated_time_in_hours 9999
```
Kill execution after N hours

```
output_file outputS.out
```
Print scattering matrix elements per scattering angle

```
pmatrix_out pmatrix.out
```
Print scattering matrix that is readable by the SIRIS (so you can use the output as an input for next round...)

```
details_out details.out
```
Print other details about the finished work

```
I_cutoff_limit 0.0000001
```
Cut off limit. When intensity of the ray goes below this limit, kill it

```
seed 0
```
Generate random seed for the PRNG (0), If nonzero, the given number will be used as a seed.

```
wavelen 6.283185307179586
```
Wavelength. Unit does not matter as long as you keep it consistent with other length parameters

```
mesh_scale 2000
```
Mesh file is scaled with this value. So, if the wavelenght is 6 nm, and if the mesh file has a sphere with radius 1 (dimensionless), SIRIS will compute 2000-nm-sized spehre with mesh_scale 2000.

```
mesh sphere.off
```
Relative path to the shape model (see below)

```
force_interaction 1
```
This is related to the diffuse scattering. Do we force every ray to interact with the diffusely scattering media, or do we let them go through. Without this the observer can see huge spike at the forward scattering direction.

```
beam_radius 250
```
The radius of the incident beam (same units as above). Negative beam size means that the entire medium is covered by the beam.

```
material1 1.0 0.0 1 1.0 0.9460290562711914 outputS_ave.out 1 ~/dists/cdfconstant04dist.txt
```
Define materialX, where the X is the number of the material. Supports 255 materials. It is important to note that material 0 is reserved for the surrounding media.  The format is
```
materialX REF_REAL REF_IMAG DIFFUSE_ON ALBEDO MEAN_FREE_PATH PATH_TO_PHASE_MATRIX EXPERIMENTAL_MFP_ON file_path
```
where REF_REAL is the refractive real part and REF_IMAG is the imaginary part. Material can have diffuse inclusions that are enabled by using 1 for DIFFUSE_ON. ALBEDO, MEAN_FREE_PATH and PATH_TO_PHASE_MATRIX are for the diffuse scatterers, whereas EXPERIMENTAL_MFP_ON is about the experimental mean free paths (see SIRIS2019 or SIRIS2020 paper). 

```
media REF_REAL REF_IMAG
```
For the surrounding medium, REF_REAL is the refractive real part and REF_IMAG is the imaginary part.

**Off file**

First, the mesh needs to be watertight and only have triangle faces. Only Off-file format is supported, and the materials are coded into the colors of the faces: Next, see an example of an OFF-file (read more from Wikipedia):   
```
OFF
# Example
 
8 6 12
 1.0  0.0 1.4142
 0.0  1.0 1.4142
-1.0  0.0 1.4142
 0.0 -1.0 1.4142
 1.0  0.0 0.0
 0.0  1.0 0.0
-1.0  0.0 0.0
 0.0 -1.0 0.0
3  0 1 2  1 0 0 
3  7 4 0  0 1 0 
3  4 5 1  0 0 0 
3  5 6 2  2 1 0 
3  3 2 6  0 0 0
3  6 5 4  1 0 0
```

Anyway, the colors are coded at lines with face information
```
3  0 1 2  1 0 0 
```
were the first '3' means that this face has 3 vertices (so triangle), then the face is made of vertices [0,1,2], which is followed by the color [1,0,0]. This color has the material coding, and this line means that the face is a surface between material 0 and 1, whereas the line
```
3  5 6 2  2 1 0 
```
has material surfaces 2 and 1. The order does not matter. So if the triangles intersected during the ray tracing have color codings [[1,0],[1,2],[1,2]], the material the ray is in is 1, because the ray tracing always starts from the material 0. Next, we encounter a triangle with encoding [1,0], meaning that the refracted ray has material coding '1' and reflected 0. The next encounter is with [1,2], and because we are in the material 1, the refracted ray must be at material 2 and so on. If the SIRIS encounters a material surface that is impossible (e.g. [[1,0],[1,2],[4,0]], where the ray encounters surface between materials 4 and 0 although the ray was in material 2, the program will stop execution and prints error about it (I hope it does). Probable causes: the object in the mesh is not watertight or has bad material coding.

### GS Gaussian spheres

Produces random Gaussian sphere shapes and discretized mesh representations for them. Can write output in Matlab, IDL IDF, Paraview VTK, and OFF.

If pubishing or distributing results that use this code, please reference to \[1\].

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

## References

1. Muinonen K, Nousiainen T, Fast P, Lumme K, and Peltoniemi JI (1996). Light scattering by Gaussian random particles: Ray optics approximation. Journal of Quantitative Spectroscopy & Radiative Transfer 55(5), 577–601. DOI:[10.1016/0022-4073(96)00003-9](https://doi.org/10.1016/0022-4073(96)00003-9).
2. Muinonen K, Nousiainen T, Lindqvist H, Muñoz O, and Videen G (2009). Light scattering by Gaussian particles with internal inclusions and roughened surfaces using ray optics. Journal of Quantitative Spectroscopy & Radiative Transfer 110, 1628–1639. DOI:[10.1016/j.jqsrt.2009.03.012](https://doi.org/10.1016/j.jqsrt.2009.03.012).
3. Lindqvist H, Martikainen J, Räbinä J, Penttilä A, and Muinonen K (2018). Ray optics in absorbing media with application to ice crystals at near-infrared wavelengths. Journal of Quantitative Spectroscopy & Radiative Transfer 217, 329–337. DOI:[10.1016/j.jqsrt.2018.06.005](https://doi.org/10.1016/j.jqsrt.2018.06.005).
4. Väisänen T, Martikainen J, and Muinonen K (2019). Scattering of light by dense particulate media in the geometric optics regime. Journal of Quantitative Spectroscopy & Radiative Transfer, in press. DOI:[10.1016/j.jqsrt.2019.106719](https://doi.org/10.1016/j.jqsrt.2019.106719).
