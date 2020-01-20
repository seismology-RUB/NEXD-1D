#-----------------------------------------------------------------------
#   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
#   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
#   Copyright 2014-2020 Marc S. Boxberg (RWTH Aachen University, GER)
#
#   This file is part of NEXD 1D.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with NEXD 2D. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------

#  General definitions
#
# Check Operating system:
OS := $(shell uname)
#
#  set the SHELL
#   MSB MSB: Actually nothing should be specified here! (or SHELL = /bin/sh if really necessary for any reason...)
SHELL = /bin/tcsh
#  Paths
bindir = ./bin
obsdir = ./obj
moduledir = ./mod
srcdir = ./src
AR = ar

# Choose your compiler:
gf = gfortran
#gf = ifort
F95 = mpif90

ifeq ($(notdir $(F95)),g95)
	FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else ifeq ($(gf),gfortran)
	# ==================================================================
	# GFORTRAN COMPILER
	# ==================================================================
	# Choose a suitable set of compiler flags for gfortran. The default 
	# flags are explained in the following...
	# required:
	#   -J$(moduledir)
	#   -fimplicit-none
	#   -ffixed-line-length-none
	#   -ffree-line-length-none
	# optional:
	#   -O3                : The compiler tries to optimize the code to 
	#                        make it faster. Level 3 is the most 
	#                        optimization available.
	#   -Og                : Enables optimizations that do not interfere
	#                        with debugging. It should be the 
	#                        optimization level of choice for the 
	#                        standard edit-compile-debug cycle.
	#   -funroll-all-loops : This might increase the speed of the code 
	#                        as well.
	#   -fbounds-check     : Add a check that the array index is within 
	#                        the bounds of the array every time an array
	#                        element is accessed.
	#                        This substantially slows down a program 
	#                        using it, but is a very useful way to find 
	#                        bugs related to arrays.
	#   -fbacktrace        : If the program crashes, a backtrace will be
	#                        produced if possible.
	#   -Wall              : gfortran will generate warnings about many 
	#                        common sources of bugs.
	#   -Wpedantic         : Generate warnings about language features 
	#                        that are supported by gfortran but are not 
	#                        part of the official Fortran 95 standard. 
	#   -Wextra            : This enables some extra warning flags that 
	#                        are not enabled by -Wall. Note, 
	#                        -Wunused-parameter (part of -Wextra) has a 
	#                        problem to correctly work with constants.h.
	# recommended flags for developers:
#	FFLAGS = -Og -J$(moduledir) -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -fbounds-check -fbacktrace -Wall -Wpedantic -Wextra
	# recommended flags for profiling (CP CP):
#	FFLAGS = -O3 -J$(moduledir) -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -funroll-all-loops -pg -g -no-pie
	# recommended flags for users:
	FFLAGS = -O3 -J$(moduledir) -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -funroll-all-loops
else ifeq ($(gf),ifort)
	# ==================================================================
	# INTEL COMPILER (ifort)
	# ==================================================================
	# Choose a suitable set of compiler flags for ifort. The default 
	# flags are explained in the following...
	# required:
	#   -module $(moduledir)
	#   -implicitnone
	#   -132
	# optional:
	#   -O3                : The compiler tries to optimize the code to 
	#                        make it faster. Level 3 is the most 
	#                        optimization available. Level 0 (-O0) is 
	#                        recommended for developers.
	#   -funroll-loops     : This might increase the speed of the code 
	#                        as well.
	#   -nowarn            : Suppresses all warnings.
	# developers:
#	FFLAGS = -O0 -module $(moduledir) -implicitnone -132
	# users:
	FFLAGS = -O3 -module $(moduledir) -implicitnone -132 -funroll-loops -nowarn
endif

obj_solver = \
	solverPro.o \
	constantsMod.o \
	parameterMod.o \
	matrixMod.o \
	mesh1DMod.o \
	gllMod.o \
	jacobiMod.o \
	vandermonde1DMod.o \
	localOperatorsMod.o \
	rosettaGammaMod.o \
	materials1DMod.o \
	outputMod.o \
	riemannflux1Dmod.o \
	rhsElastic1DMod.o \
	rhsPoroelastic1DMod.o \
	rhsSlipInterface1DMod.o \
	errorMessage.o \
	realloc.o \
	sourcesMod.o \
	slipInterfaceMod.o \
	genericMaterial.o \
	analyticalSolutionMod.o \
	slopeLimiterMod.o \
	receiverMod.o \
	boundaryConditionsMod.o \
	calendar.o \
	convert_time.o
	
obj_analyticalSol = \
	analyticalSolution.o \
	constantsMod.o \
	parameterMod.o \
	matrixMod.o \
	slipInterfaceMod.o \
	errorMessage.o \
	materials1DMod.o \
	realloc.o
	
obj_test = testprogramm.o \
#-------------------------------------------------------
#  Direcory search
#
vpath %.o $(obsdir)
vpath %.f90 $(srcdir) ./src/include
vpath %.f $(srcdir) ./src/include
vpath %.c $(srcdir) ./src/include

#--------------------------------------------------------
#  additional directories to be searched for module or include dependencies
#  default is search in ./ only
#
DEPDIRS = $(srcdir) ./src/include
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(gf) -c $(FFLAGS) $< -o $(obsdir)/$@
%.o: %.f
	$(gf) -c $(FFLAGS) -fimplicit-none -ffixed-line-length-132 $< -o $(obsdir)/$@
%.o: %.c	
	gcc -c $(CFLAGS) $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
obstringtest = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Library paths
ifeq ($(OS), Linux)
	pgplot = -lpgplot -L/usr/lib -lpng -lz -lX11
endif
ifeq ($(OS), Darwin)
	pgplot = -lpgplot -L/usr/local/lib -lpng -lz -L/usr/X11/lib -lX11
endif
la = -llapack -lblas
#lib = metis-4.0.3/libmetis.a 
#---------------------------------------------------------
.PHONY: all
#
#  create dependencies on modules 
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make.incdep:
	./makeDepFromUseInclude.py $(DEPDIRS) > $@
-include make.incdep
#
cleanlsi:	
	rm -f $(bindir)/solver $(moduledir)/*.mod $(obsdir)/*.o make.incdep
#
clean: 
	rm -f $(bindir)/solver $(bindir)/analytical $(moduledir)/*.mod $(obsdir)/*.o make.incdep
#
solver: $(obj_solver)
	$(gf) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(la) $(pgplot)
	
analytical: $(obj_analyticalSol)
	$(gf) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(la)
	
test: $(obj_test)
	$(gf) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(la)
	
bin:
	mkdir -p $(bindir)
#
mod:
	mkdir -p $(moduledir)
#
obj:
	mkdir -p $(obsdir)
#
required: mod obj bin
#
DEFAULT = \
	required \
	solver
#
default: $(DEFAULT)
#
all: default
