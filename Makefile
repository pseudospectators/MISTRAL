# Makefile for fsi and mhd codes. See README for necessary environment
# variables.

# If no HDF5 library is set, compile without HDF5 support
#ifndef HDF_ROOT
#HDF5FLAG = no
#endif

# HDF5 flag set to 'yes' by default
ifndef HDF5FLAG
HDF5FLAG = yes
endif

# Non-module Fortran files to be compiled:
FFILES = rhs.f90 fluid_time_step.f90 init_fields.f90 \
	mask.f90 mask_fsi.f90 time_step.f90 \
	init_fields_fsi.f90 integrals.f90 params.f90 \
	runtime_control.f90 drag.f90 \
	draw_plate.f90 draw_sphere.f90 \
  rotation_matrices.f90  basic_file_routines.f90 \
  add_channel.f90 add_cavity.f90 \
  noncircular_cylinder.f90 draw_flexible_plate.f90 implicit.f90

ifeq ($(HDF5FLAG),yes)
FFILES += save_fields.f90 postprocessing.f90
else
FFILES += save_fields_nohdf5.f90
endif

# Object and module directory:
OBJDIR=obj
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = vars.f90 diff.f90 cof_p3dfft.f90 solid_solver.f90 \
	interpolation.f90 basic_operators.f90 insects.f90 ghostpoints.f90 \
	ini_files_parser.f90 ini_files_parser_mpi.f90 helpers.f90 hdf5_wrapper.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = src
VPATH += :src/inicond:src/inicond/hyd:src/inicond/mhd:src/inicond/scalar
VPATH += :src/geometry:src/geometry/hyd:src/file_io
VPATH += :src/insects:src/solid_solver

# Set the default compiler if it's not already set, make sure it's not F77.
ifndef FC
FC = mpif90
endif
ifeq ($(FC),f77)
FC = mpif90
endif
ifeq ($(FC),sxf90)
FC = sxmpif90
endif

#SX compiler
ifeq ($(FC),sxmpif90)
FFLAGS += -I$(OBJDIR)
FFLAGS += -R2 -Wf"-pvctl fullmsg" -ftrace -f2003
PPFLAG =
else

# GNU compiler
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
FFLAGS += -Wall # warn for unused and uninitialzied variables
#FFLAGS += -Werror # warnings are errors
FFLAGS += -pedantic
PPFLAG = -cpp #preprocessor flag
# Debug flags for gfortran:
#FFLAGS += -Wuninitialized -O -fimplicit-none -fbounds-check -g -ggdb
endif

# Intel compiler
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
PPFLAG = -fpp #preprocessor flag
DIFORT = -DIFORT # define the IFORT variable
FFLAGS += -module $(OBJDIR) # specify directory for modules.
FFLAGS += -vec_report0
#FFLAGS += -mcmodel=medium -i-dynamic
# debug flags for ifort:
#FFLAGS += -g -pg
#FFLAGS +=-CB -traceback -debug extended -fpe0 -warn -check all -gen-interfaces -warn interfaces
endif

#IBM compiler
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
FFLAGS += -qmoddir=$(OBJDIR)
FFLAGS += -I$(OBJDIR)
DIFORT =-WF,-DTURING # this defines the TURING with the IBM compiler
PPFLAG =-qsuffix=cpp=f90  #preprocessor flag
endif

endif

PROGRAMS = mistral

# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include

# Linker flags
ifeq ($(FC),sxmpif90)
LDFLAGS += -lblas
endif

LDFLAGS =
ifeq ($(HDF5FLAG),yes)
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB)
LDFLAGS += -lhdf5_fortran -lhdf5
endif
LDFLAGS += -llapack
ifeq ($(FC),sxmpif90)
LDFLAGS += -lblas
else
LDFLAGS += -lz -ldl
endif
LDFLAGS += -lm

ifeq ($(HDF5FLAG),yes)
FFLAGS += -I$(HDF_INC)
endif
FFLAGS += $(PPFLAG) $(DIFORT)

# Both programs are compiled by default.
all: directories $(PROGRAMS)

# Compile main programs, with dependencies.
mistral: mistral.f90 $(MOBJS) $(OBJS)
	rm -rf mistral
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).
$(OBJDIR)/vars.o: vars.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/diff.o: diff.f90 $(OBJDIR)/vars.o $(OBJDIR)/ghostpoints.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/cof_p3dfft.o: cof_p3dfft.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/insects.o: insects.f90 $(OBJDIR)/vars.o \
	body_geometry.f90 body_motion.f90 rigid_solid_time_stepper.f90 wings_geometry.f90 wings_motion.f90 stroke_plane.f90 \
  kineloader.f90 $(OBJDIR)/helpers.o $(OBJDIR)/ini_files_parser_mpi.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/solid_solver.o: solid_solver.f90 $(OBJDIR)/vars.o  $(OBJDIR)/interpolation.o $(OBJDIR)/basic_operators.o $(OBJDIR)/insects.o \
	mouvement.f90 integrate_position.f90 init_beam.f90 save_beam.f90 BeamForces.f90 plate_geometry.f90 $(OBJDIR)/ghostpoints.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/interpolation.o: interpolation.f90 $(OBJDIR)/vars.o $(OBJDIR)/basic_operators.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/basic_operators.o: basic_operators.f90 $(OBJDIR)/vars.o $(OBJDIR)/cof_p3dfft.o $(OBJDIR)/diff.o \
	$(OBJDIR)/ghostpoints.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ghostpoints.o: ghostpoints.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser_mpi.o: ini_files_parser_mpi.f90 $(OBJDIR)/vars.o $(OBJDIR)/ini_files_parser.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser.o: ini_files_parser.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/helpers.o: helpers.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/hdf5_wrapper.o: hdf5_wrapper.f90 $(OBJDIR)/vars.o $(OBJDIR)/helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR)/*.o $(OBJDIR)/*.mod *.L *.mod mistral.o a.out

tidy:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod *.L *.mod mistral.o a.out

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}
	@echo 'Build with HDF5 support:' $(HDF5FLAG)

${OBJDIR}:
	mkdir -p ${OBJDIR}
