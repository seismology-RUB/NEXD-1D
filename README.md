# NEXD 1D

NEXD is a software package for high order simulation of seismic waves using the nodal discontinuous Galerkin method.
NEXD 1D is the 1D solver of this software package.

NEXD uses the Nodal Discontinuous Galerkin Method (NDG) to calculate synthetic seismograms and supports a high order spatial approximation of the wave field. The simulation tool NEXD is developed with the capability of simulating elastic and anelastic wave fields for seismic experiments for one-, two- and three- dimensional settings. The implementation of poroelasticity and simulation of slip interfaces as representation of fractures are currently in progress and are working for the one- and two-dimensional part.

## Authors and License

NEXD 1D and some of its components, as well as documentation and some examples
are available under terms of the [GNU General Public License](LICENSE) (version 3 or higher)
on [github](https://github.com/seismology-RUB/NEXD-1D).
Please find contact addresses [there](https://github.com/seismology-RUB), or visit 
http://www.rub.de/nexd in case you want to get in touch with the authors. If you 
encounter any problems installing or using the software, it will be helpful to 
open (or add to) an "issues" topic at the [github repository](https://github.com/seismology-RUB/NEXD-1D).

The main authors are Thomas Möller (Ruhr-Universität Bochum, Germany) and Marc S. Boxberg (RWTH Aachen University, Germany).


## Documentation

A user manual is in preparation and will be available with the next release.


## Examples

NEXD 1D comes with three examples: an example for elastic wave propagation (simulations/example), an example for poroelastic wave propagation (simulations/example_poroelastic) and an example for including fractures (simulations/example_fracture). Please use these examples to familiarize yourself with the software.


## Requirements

* GNU Make
* Fortran compiler (sufficient standard)
* LAPACK libraries
* pgplot [http://www.astro.caltech.edu/~tjp/pgplot](http://www.astro.caltech.edu/~tjp/pgplot)


## Installation

0. If not yet done, you should download the source code of the NEXD 1D main package by cloning the master branch of the NEXD 1D repository on gitHub.com:
     ```
     git clone --depth 1 --branch master https://github.com/seismology-RUB/NEXD-1D
     ```

1. Install all software dependencies.

2. Adjust the software to your system and personal requirements by changing the [Makefile](Makefile) appropriately (e.g., set your compiler).

3. Run command
     ```
     make all
     ```
   from installation path `NEXD-1D/` to compile NEXD 1D.
