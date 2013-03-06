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
  timingmod.o

LIBRARIES = GGRID/ggrid.a MISC/misc.a

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
version.inc: version.py
	@python2 -B version.inc.py

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
banner.o: version.inc
dealloc_all.o: bbxsmod.o gasgridmod.o ionsmod.o mpimod.o
read_bbxs_data.o: bbxsmod.o ionsmod.o physconstmod.o timingmod.o
read_restart_file.o: inputparmod.o miscmod.o gasgridmod.o tempcorrmod.o timestepmod.o

supernu.o: bfxsmod.o ffxsmod.o gasgridmod.o gstructmod.o inputparmod.o ionsmod.o mpimod.o physconstmod.o timestepmod.o timingmod.o

#
#-- LIBRARIES
$(SUBDIRS):
	$(MAKE) -C $@
#
#-- PROGRAMS
supernu: $(MODULES) snipper.o banner.o \
  read_bbxs_data.o read_restart_file.o write_restart_file.o \
  dealloc_all.o $(LIBRARIES)
