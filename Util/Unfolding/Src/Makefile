# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

# Makefile for unfold

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

.PHONY: default what dep version  clean

override WITH_NCDF=
override WITH_FLOOK=

default: what unfold 

VPATH:=$(shell pwd)/../../../Src
OBJDIR=Obj
MAIN_OBJDIR:=$(shell cd -P -- ../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)

# Set the default non-mpi routines
FC_DEFAULT := $(FC)
FC_SERIAL  ?= $(FC_DEFAULT)

# Filter out non-needed libraries
COMP_LIBS := $(filter-out libncdf.a, $(COMP_LIBS))
COMP_LIBS := $(filter-out libfdict.a, $(COMP_LIBS))

# For removing any -ipo compilation in sub directories...
# Later compilers are also having this enabled, and it lets you
# specify what to not be used in libraries.
# This is so because inter-procedural compilation should NOT be performed in
# libraries (it removes names when a routine can be fully moved into others)
IPO_FLAG?= -ipo

SIESTA_SYS=$(SIESTA_ARCH)

what:
	@echo
	@echo "Compilation architecture to be used: $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo


NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src
INCFLAGS:=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSML_INCFLAGS) $(PSOP_INCFLAGS) \
          $(INCFLAGS)

XC = $(GRIDXC_LIBS)
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a

INCFLAGS += $(GRIDXC_INCFLAGS)

SYSOBJ=$(SYS).o

SIESTA_OBJS= precision.o m_io.o alloc.o memory_log.o memory.o pxf.o \
             moreParallelSubs.o m_mpi_utils.o atmfuncs.o atm_types.o \
             radial.o reclat.o idiag.o minvec.o \
             parallel.o io.o files.o cellsubs.o sorting.o \
	     timestamp.o m_wallclock.o \
             diag.o diag_option.o bloch_unfold.o \
	     debugmpi.o find_kgrid.o units.o posix_calls.o xml.o \
             interpolation.o spher_harm.o basis_types.o atmparams.o basis_io.o \
             atom_options.o chemical.o \
             get_kpoints_scale.o siesta_geom.o \
             radfft.o bessph.o m_fft_gpfa.o \
             broadcast_basis.o dftu_specs.o basis_specs.o \
             atom.o hamann.o m_filter.o \
             m_cite.o \
             coor.o chkdim.o \
             zmatrix.o m_cell.o  periodic_table.o volcel.o digcel.o \
             arw.o redcel.o sys.o parallelsubs.o spatial.o \
             domain_decom.o printmatrix.o schecomm.o \
             pspltm1.o mmio.o qsort.o \
             sparse_matrices.o

SIESTA_OBJS += m_walltime.o timer_tree.o m_timer.o timer.o
SIESTA_OBJS += m_spin.o m_vee_integrals.o


# Add specific class objects
# class Objects
SIESTA_OBJS += class_OrbitalDistribution.o \
	class_Sparsity.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        class_Fstack_Pair_Geometry_SpData2D.o \
        class_Pair_Geometry_SpData2D.o \
        class_Geometry.o \
        m_uuid.o object_debug.o

# As the classes are based on inputs we force the dependency to
# alloc
class_Sparsity.o class_Data1D.o class_Data2D.o: alloc.o

DEP_OBJS = $(SIESTA_OBJS) $(TS_OBJS) $(SYSOBJ) $(EXTRA)

# To separate unfold objects from siesta/transiesta
UNFOLD_OBJS  = array.o hsx_m.o
UNFOLD_OBJS += unfold.o 

# Handlers and resolvers of last resort
UNFOLD_OBJS += handlers.o
#
# Extract information of the compiler version used.
# This is very basic but for intel, gcc, pgi and open64
# this at least contain the version string of the compiler.
# Possibly later the used compiler version may be extracted.
COMPILER_VERSION := $(shell $(FC) --version | sed '/^[ ]*$$/d;q')
ifeq ($(COMPILER_VERSION),)
  # Just try and capture the case were it can not figure out anything
  COMPILER_VERSION := Unknown version
endif

# Interfaces to libraries
MPI_MAKEFILE=$(VPATH)/MPI/Makefile
libmpi_f90.a: 
	-mkdir -p MPI
	(cd MPI ; \
	$(MAKE) -f $(MPI_MAKEFILE) -j 1 "VPATH=$(VPATH)/MPI" \
	"MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=$(ARCH_MAKE)" \
	"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module_built)

FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf -I../MPI $(INCFLAGS)
$(FDF): $(MPI_INTERFACE)
	-mkdir -p fdf
	(cd fdf ; \
	$(MAKE) -f $(FDF_MAKEFILE) -j 1 "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
	"ARCH_MAKE=$(ARCH_MAKE)" \
	"INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module)

#######################
#
# Linear algebra routines
#
libsiestaBLAS.a: $(VPATH)/Libs/blas.F
	@echo "==> Compiling libsiestaBLAS.a in Libs..."
	-mkdir -p Libs
	@(cd Libs ; $(MAKE) -f $(VPATH)/Libs/Makefile -j 1 \
		"VPATH=$(VPATH)/Libs" \
		"ARCH_MAKE=$(ARCH_MAKE)" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaBLAS.a)
	@cp Libs/libsiestaBLAS.a .
libsiestaLAPACK.a: $(VPATH)/Libs/lapack.F
	@echo "==> Compiling libsiestaLAPACK.a in Libs..."
	-mkdir -p Libs
	@(cd Libs ; $(MAKE) -f $(VPATH)/Libs/Makefile -j 1 \
		"VPATH=$(VPATH)/Libs" \
		"ARCH_MAKE=$(ARCH_MAKE)" \
		"FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" libsiestaLAPACK.a)
	@cp Libs/libsiestaLAPACK.a .


# All objects require the FDF library
# In practice they do not but it greatly simplifies the logic
# here.
# This will setup the correct compilation sequence as
# FDF => MPI
$(DEP_OBJS): $(FDF) $(COMP_LIBS) 

unfold: $(MPI_INTERFACE) $(FDF) $(COMP_LIBS) $(DEP_OBJS) $(UNFOLD_OBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o unfold \
              $(DEP_OBJS) $(UNFOLD_OBJS) $(MPI_INTERFACE) \
              $(FDF) $(XC) $(NCPS) $(PSOP) $(PSML_LIBS) $(XMLF90_LIBS) $(COMP_LIBS) $(LIBS)

.PHONY: clean
clean: 
	@echo "==> Cleaning object, library, and executable files"
	-rm -f unfold *.o *.a *.mod
	-rm -rf ./fdf
	-rm -rf ./MPI
	-rm -rf ./Libs

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(VPATH)/class*.T90 \
		$(UNFOLD_OBJS:.o=.f90) $(UNFOLD_OBJS:.o=.F90) \
		$(UNFOLD_OBJS:.o=.f) $(UNFOLD_OBJS:.o=.F) || true

	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# Assert that FDF is compiled on before hand
$(DEP_OBJS): $(FDF)

# DO NOT DELETE THIS LINE - used by make depend
arw.o: alloc.o parallel.o precision.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atm_types.o atmparams.o atom_options.o basis_specs.o basis_types.o
atom.o: hamann.o interpolation.o m_filter.o periodic_table.o precision.o
atom.o: radial.o sys.o
atom_options.o: sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o radial.o sys.o xml.o
basis_specs.o: atom_options.o basis_types.o chemical.o periodic_table.o
basis_specs.o: precision.o sys.o
basis_types.o: alloc.o atmparams.o precision.o sys.o
bessph.o: precision.o sys.o
bloch_unfold.o: units.o
broadcast_basis.o: atm_types.o chemical.o dftu_specs.o m_mpi_utils.o m_spin.o
broadcast_basis.o: parallel.o radial.o
cellsubs.o: precision.o
chemical.o: parallel.o precision.o sys.o
chkdim.o: sys.o
class_Data1D.o: alloc.o
class_Data2D.o: alloc.o
class_Data3D.o: alloc.o
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Geometry.o: alloc.o
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_SpData1D.o: class_Data1D.o
class_SpData1D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData2D.o: class_Data2D.o
class_SpData2D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData3D.o: class_OrbitalDistribution.o class_Sparsity.o
class_Sparsity.o: alloc.o
class_TriMat.o: alloc.o
coor.o: alloc.o atmfuncs.o parallel.o precision.o siesta_geom.o sys.o units.o
coor.o: zmatrix.o
debugmpi.o: parallel.o
dftu_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
dftu_specs.o: basis_specs.o basis_types.o interpolation.o m_cite.o m_spin.o
dftu_specs.o: m_vee_integrals.o parallel.o precision.o radial.o sys.o units.o
diag.o: alloc.o diag_option.o parallel.o precision.o sys.o
diag_option.o: parallel.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
get_kpoints_scale.o: parallel.o siesta_geom.o units.o
idiag.o: parallel.o sys.o
io.o: m_io.o
m_cell.o: precision.o siesta_geom.o units.o
m_filter.o: bessph.o precision.o radfft.o sys.o
m_io.o: sys.o
m_mpi_utils.o: precision.o sys.o
m_spin.o: precision.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o sys.o
m_vee_integrals.o: precision.o sys.o
m_wallclock.o: m_walltime.o
memory.o: memory_log.o
memory_log.o: m_io.o parallel.o precision.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o sys.o
parallelsubs.o: domain_decom.o parallel.o precision.o spatial.o sys.o
printmatrix.o: alloc.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
redcel.o: parallel.o precision.o sys.o
schecomm.o: alloc.o
siesta_geom.o: precision.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
sparse_matrices.o: precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o sys.o
timer.o: m_timer.o parallel.o sys.o timer_tree.o
timer_tree.o: m_walltime.o
units.o: precision.o
xml.o: precision.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o sys.o units.o
array.o: precision.o sys.o
unfold.o: alloc.o array.o atmfuncs.o basis_io.o basis_types.o cellsubs.o
unfold.o: diag_option.o get_kpoints_scale.o hsx_m.o m_io.o m_mpi_utils.o
unfold.o: m_timer.o memory_log.o parallel.o parallelsubs.o precision.o radfft.o
unfold.o: siesta_geom.o spher_harm.o sys.o
aux_proj.o: atom.o
bloch_unfold_m.o: bloch_unfold.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_fstack_pair_geometry_dspdata2d.o: class_Fstack_Pair_Geometry_SpData2D.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_ldata1d.o: class_Data1D.o
class_lspdata1d.o: class_SpData1D.o
class_pair_geometry_dspdata2d.o: class_Pair_Geometry_SpData2D.o
class_sdata1d.o: class_Data1D.o
class_sspdata1d.o: class_SpData1D.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
m_array.o: array.o
m_bessph.o: bessph.o
m_diag.o: diag.o
m_diag_option.o: diag_option.o
m_find_kgrid.o: find_kgrid.o
m_get_kpoints_scale.o: get_kpoints_scale.o
m_hamann.o: hamann.o
m_minvec.o: minvec.o
m_object_debug.o: object_debug.o
m_radfft.o: radfft.o
m_timer_tree.o: timer_tree.o
m_timestamp.o: timestamp.o
mtprng.o: m_uuid.o
printmat.o: printmatrix.o
t_spin.o: m_spin.o
timer_options.o: timer.o
