# Makefile for kbgen
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: clean
#
default: psml2psf
#
MAIN_OBJDIR=.
TOPDIR=.

VPATH:=$(TOPDIR)/Src:$(TOPDIR)/Pseudo/converters/psml2psf

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#

INCFLAGS+= $(PSML_INCFLAGS) $(NCPS_INCFLAGS)
INCFLAGS+= $(XMLF90_INCFLAGS) $(GRIDXC_INCFLAGS)
INCFLAGS+= $(LIBSYS_INCFLAGS)
#
OBJS = 	psml2psf.o 
#
#
MOD_OBJS= m_getopts.o
#
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): 
#
what:
	echo "VPATH: " $(VPATH)
psml2psf:   $(NCPS)  $(ALL_OBJS) $(LIBSYS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psml2psf $(ALL_OBJS) \
        $(NCPS) $(PSML_LIBS) $(XMLF90_LIBS) $(LIBSYS)
#
psml2psf.o: m_getopts.o 
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f psml2psf *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
PROGS:= psml2psf
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin



