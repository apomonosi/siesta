#
# Defs for a GCC system, with all library symbols loaded from modules
#
# ml flook netcdf openmpi
# ml scalapack      # defines SCALAPACK_LIBS
# ml lapack   (Mac with homebrew: ml veclibfort)  # defines LAPACK_LIBS
# ml fftw
# ml elpa
# ml gridxc-multi libpsml
#--------------------------------------------------------
# Use these symbols to request particular features
# To turn on, set '=1'.
#--------------
#
WITH_MPI=1
WITH_ELPA=1
WITH_FLOOK=1
WITH_NETCDF=1
WITH_NCDF=1
WITH_LEGACY_GRIDXC_INSTALL=
WITH_GRID_SP=

#-------------
# Define compiler names and flags
#
FC_PARALLEL=mpif90
FC_SERIAL=gfortran
#
FPP = $(FC_SERIAL) -E -P -x c
FFLAGS= -O2 -g 
FFLAGS_DEBUG= -g -O0 #-fcheck=all
RANLIB=echo
# ----------------------------------------------------------
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include $(SELF_DIR)build.mk
