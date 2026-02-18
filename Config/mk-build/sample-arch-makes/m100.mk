# Arch.make file for Marconi-100 at CINECA
# with the gcc suite, including GPU support
#
#-------------------------------------------
# Example of modules to be loaded:
#
# load("gnu/8.4.0")
# load("spectrum_mpi/10.3.1--binary")
# load("blas/3.8.0--gnu--8.4.0")
# load("lapack/3.9.0--gnu--8.4.0")
# load("scalapack/2.1.0--spectrum_mpi--10.3.1--binary")
# load("essl/6.2.1--binary")
# prepend_path("LD_LIBRARY_PATH","/cineca/prod/opt/compilers/xl/16.1.1/binary/lib")
# load("cuda/11.0")
# conflict("DEVELOP")
# LmodMessage("display GCC devel environment")
# setenv("MPI_ROOT","/cineca/prod/opt/compilers/spectrum_mpi/10.3.1/binary")
# setenv("SCALAPACK_LIBS","-L/cineca/prod/opt/libraries/scalapack/2.1.0/spectrum_mpi--10.3.1--binary/lib -lscalapack")
# setenv("LAPACK_LIBS","-L/cineca/prod/opt/libraries/lapack/3.9.0/gnu--8.4.0/lib -llapack -L/cineca/prod/opt/libraries/essl/6.2.1/binary/lib64 -lessl")
# setenv("CUDA_ROOT","/cineca/prod/opt/compilers/cuda/11.0/none/lib64")
# setenv("CUDA_LIBS","-L/cineca/prod/opt/compilers/cuda/11.0/none/lib64 -lcublas -lcudart")
# prepend_path("LD_LIBRARY_PATH","/cineca/prod/opt/compilers/cuda/11.0/none/lib64")
# setenv("OMPI_FC","gfortran")

# (+ appropriate libpsml, gridxc-multi, elpa, etc)
#

WITH_ELPA=1
WITH_FLOOK=
WITH_MPI=1
WITH_NETCDF=
WITH_NCDF=
WITH_LEGACY_GRIDXC_INSTALL=
WITH_GRID_SP=
#

#  ml load gcc_env elsi libpsml xmlf90 gridxc-multi
#
#NETCDF_ROOT=$(NETCDF_HOME)            
#NETCDF_FORTRAN_ROOT=$(NETCDF_HOME)
#HDF5_LIBS=-L/apps/HDF5/1.8.20/GCC/OPENMPI/lib -lhdf5_hl -lhdf5 -lcurl -lz
#
#--- From gcc_env module:
#SCALAPACK_LIBS=
#LAPACK_LIBS=
#CUDA_LIBS=
#-----------------------
#FFTW_ROOT=
# Needed for PEXSI (ELSI) support
#LIBS_CPLUS=-lstdc++ -lmpi_cxx $(CUDA_LIBS)    # Some linux systems with OpenMPI
#LIBS_CPLUS=-lstdc++ $(CUDA_LIBS)
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC_PARALLEL=mpif90
FC_SERIAL=gfortran
FPP = $(FC_SERIAL) -E -P -x c
FFLAGS = -O2 
FFLAGS_DEBUG= -g -O0
RANLIB=echo

#
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include $(SELF_DIR)build.mk
