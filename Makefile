include Makefile.compiler
include Makefile.inc

########################################################################
# MACROS
########################################################################
PROGRAMS = supernu

MODULES = kindmod.o \
  elemdatamod.o physconstmod.o mpimod.o \
  miscmod.o \
  ionsmod.o ffxsmod.o bfxsmod.o bbxsmod.o \
  inputparmod.o timestepmod.o \
  gasgridmod.o \
  particlemod.o \
  timingmod.o

OBJFILES = grids.o initials.o xsections.o sourcenumbers.o vacancies.o interior_source.o advance.o material_update.o write_output.o diffusion1.o transport1.o dealloc_all.o

LIBRARIES = GASGRID/gasgrid.a MISC/misc.a
LIBRARIES = MISC/misc.a

SUBDIRS = $(dir $(LIBRARIES))

########################################################################
# TARGETS
########################################################################
# Utility targets (ignore corresponding file names)
.PHONY: all put-intrepid clean cleandirs cleanall run $(SUBDIRS)

all: $(MODULES) $(SUBDIRS) $(PROGRAMS)

clean: cleandirs
	rm -f *.o *.a *.mod *.MOD version.inc
	rm -f $(PROGRAMS)

cleandirs:
	for d in $(SUBDIRS); do ($(MAKE) -C $$d clean); done

run:
	BINDIR=`pwd`
	mkdir ../run || rm -rf ../run/*
	cd ../run; ln -s ${CURDIR}/Data/* .
	cd ../run; ln -s ${CURDIR}/Input/* .
	cd ../run; ln -s ${CURDIR}/supernu .
	cd ../run; supernu

########################################################################
# EXPLICIT RULES
########################################################################
#-- AUTOMATICALLY GENERATED INC FILES
version.inc: version.py
	@python2 -B version.inc.py

#
#-- MODULES
bbxsmod.o: elemdatamod.o miscmod.o physconstmod.o
gasgridmod.o: kindmod.o inputparmod.o
inputparmod.o: kindmod.o miscmod.o
ionsmod.o: miscmod.o physconstmod.o
miscmod.o: MISC/warn.f MISC/lcase.f
mpimod.o: gasgridmod.o inputparmod.o timestepmod.o timingmod.o
timestepmod.o: inputparmod.o physconstmod.o
physconstmod.o: kindmod.o
particlemod.o: kindmod.o

#
#-- OBJ FILES
#-- note: prerequisites don't need to include modules as these are always built first
banner.o: version.inc

driver.o: $(OBJFILES)

#
#-- LIBRARIES
$(SUBDIRS):
	$(MAKE) -C $@
#
#-- PROGRAMS
supernu: $(MODULES) $(OBJFILES) driver.o banner.o $(LIBRARIES)
