include Makefile.compiler
include Makefile.inc

########################################################################
# MACROS
########################################################################
PROGRAMS = supernu

MODULES = \
  elemdatamod.o physconstmod.o mpimod.o \
  miscmod.o \
  ionsmod.o ffxsmod.o bfxsmod.o bbxsmod.o \
  inputparmod.o timestepmod.o \
  gasgridmod.o \
  particlemod.o \
  timingmod.o

OBJFILES = sourcenumbers.o vacancies.o interior_source.o advance.o \
 write_output.o diffusion1.o transport1.o \
 read_bbxs_data.o dealloc_all.o

LIBRARIES = GASGRID/gasgrid.a MISC/misc.a
SUBDIRS = $(dir $(LIBRARIES))

TESTS = Testsuite/simple1 Testsuit/simple2

VERSIONPY = $(wildcard version.py)

########################################################################
#-- CONSTANTS
########################################################################
RUNDIR = $(CURDIR)/Run

########################################################################
# TARGETS
########################################################################
# Utility targets (ignore corresponding file names)
.PHONY: all put-intrepid clean cleandirs cleanall run tar testsuite $(SUBDIRS)

all: $(MODULES) $(SUBDIRS) $(PROGRAMS)

clean: cleandirs
	rm -f *.o *.a *.mod *.MOD version.inc
	rm -f $(PROGRAMS)

cleandirs:
	for d in $(SUBDIRS); do ($(MAKE) -C $$d clean); done

run: all
	mkdir $(RUNDIR) || rm -f $(RUNDIR)/*
	cd $(RUNDIR); ln -s $(CURDIR)/Data/* .
	cd $(RUNDIR); ln -s $(CURDIR)/Input/* .
	cd $(RUNDIR); ln -s $(CURDIR)/supernu .
	cd $(RUNDIR); ./supernu

testsuite: $(TESTS)
	#stub

########################################################################
# EXPLICIT RULES
########################################################################
#-- AUTOMATICALLY GENERATED INC FILES
ifeq ($(VERSIONPY), version.py)
  version.inc: $(VERSIONPY)
	@python2 -B version.inc.py
else
  version.inc: version_dummy.inc
	@cp -vu version_dummy.inc version.inc
endif

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
#-- note: prerequisites don't need to include modules as these are always built first
banner.o: version.inc

supernu.o: $(OBJFILES)

#
#-- LIBRARIES
$(SUBDIRS):
	$(MAKE) -C $@
#
#-- PROGRAMS
supernu: $(MODULES) $(OBJFILES) supernu.o banner.o $(LIBRARIES)
