# Makefile for kbgen
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: clean
#
default: psop
#
MAIN_OBJDIR=.
TOPDIR=.

VPATH:=$(TOPDIR)/Src:$(TOPDIR)/Pseudo/vnl-operator

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#

INCFLAGS += $(PSML_INCFLAGS) $(NCPS_INCFLAGS) $(PSOP_INCFLAGS)
INCFLAGS += $(XMLF90_INCFLAGS) $(GRIDXC_INCFLAGS)
INCFLAGS += $(LIBSYS_INCFLAGS) 
#
OBJS = 	psop_options.o\
	psop.o m_kb.o  \
        periodic_table.o \
        check_grid.o m_uuid.o \
        m_libxc_list.o m_libxc_sxc_translation.o \
        m_siestaxc_list.o
#
MOD_OBJS= m_getopts.o
#
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
XC=$(GRIDXC_LIBS)
#
psop:   $(NCPS) $(PSOP) $(ALL_OBJS) $(LIBSYS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psop $(ALL_OBJS) \
              $(NCPS) $(PSOP) $(PSML_LIBS) $(XC) $(XMLF90_LIBS) $(LIBSYS)
#
psop.o: psop_options.o m_uuid.o m_kb.o
psop.o: m_getopts.o periodic_table.o
psop.o: m_libxc_sxc_translation.o 
m_libxc_sxc_translation.o: m_libxc_list.o  m_siestaxc_list.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f psop *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
PROGS:= psop
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin




