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

OBJFILES = grids.o initials.o xsections.o sourcenumbers.o vacancies.o interior_source.o advance.o material_update.o write_output.o globalallocations.o diffusion1.o transport1.o

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

########################################################################
# EXPLICIT RULES
########################################################################
#-- AUTOMATICALLY GENERATED INC FILES
#version.inc: version.py
#	@python2 -B version.inc.py

#
#-- MODULES
bbxsmod.o: elemdatamod.o miscmod.o physconstmod.o
gasgridmod.o: inputparmod.o
inputparmod.o: miscmod.o
ionsmod.o: miscmod.o physconstmod.o
miscmod.o: MISC/warn.f MISC/lcase.f
mpimod.o: gasgridmod.o inputparmod.o timestepmod.o timingmod.o
timestepmod.o: inputparmod.o physconstmod.o

#
#-- OBJ FILES
#banner.o: version.inc
dealloc_all.o: bbxsmod.o gasgridmod.o ionsmod.o mpimod.o
read_bbxs_data.o: bbxsmod.o ionsmod.o physconstmod.o timingmod.o

#grids: gasgridmod.o timestepmod.o imputparmod.o
#initials: gasgridmod.o particlemod.o timestepmod.o physconstmod.o inputparmod.o
#xsections: gasgridmod.o timestepmod.o physconstmod.o inputparmod.o
#sourcenumbers: gasgridmod.o timestepmod.o particlemod.o physconstmod.o inputparmod.o
#vacancies:
#interior_source:
#advance:
#material_update:
#write_output:
#globalallocations:

supernu.o: bfxsmod.o ffxsmod.o gasgridmod.o inputparmod.o ionsmod.o mpimod.o physconstmod.o timestepmod.o timingmod.o \
  $(OBJFILES)

#
#-- LIBRARIES
$(SUBDIRS):
	$(MAKE) -C $@
#
#-- PROGRAMS
#supernu: $(MODULES) supernu.o banner.o \
  read_bbxs_data.o \
  dealloc_all.o $(LIBRARIES)
#supernu: $(MODULES) $(OBJFILES) supernu.o banner.o $(LIBRARIES)
supernu: $(MODULES) $(OBJFILES) supernu.o $(LIBRARIES)