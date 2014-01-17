########################################################################
# Makefile.compiler: compiler flags (compiler specific)
# Makefile.include : general implicit make rules
# Makefile.depend  : explicit dependencies for object files.  This file
#   is generated automatically by the Makefile.depend: rule, using the
#   depend.sh bash script.  This script extracts the "^ use " and
#   "^ include " lines from the source files.  Note that only the lower
#   case use and includes are extracted.  This allows to use the upper
#   case USE and INCLUDE statements to bypass this scan.  This is used
#   in two situations:
#   1) If the 'use'd module is defined in the same source file so that
#      the corresponding .mod file is no make dependency.
#   2) If the 'INCLUDE'd file is not in the $(TOP)/Include directory,
#      for example mpi.h
########################################################################
TOP := $(CURDIR)

include Makefile.compiler
include Makefile.inc

PROGRAMS := supernu

MODULES := gammaprofmod.o \
  elemdatamod.o physconstmod.o mpimod.o \
  miscmod.o \
  ionsmod.o ffxsmod.o bfxsmod.o bbxsmod.o \
  inputparmod.o timestepmod.o \
  gasgridmod.o inputstrmod.o \
  particlemod.o \
  timingmod.o manufacmod.o

FILES := sourcenumbers.o vacancies.o boundary_source.o interior_source.o particle_advance.o \
  write_output.o diffusion1.o transport1.o read_wlgrid.o \
  read_bbxs_data.o restart_file.o dealloc_all.o specint.o initialnumbers.o binsrch.o \
  initialnumbers.o initial_particles.o energy_check.o tau_update.o \
  banner.o

LIBRARIES := GASGRID/gasgrid.a MISC/misc.a
SUBDIRS := $(dir $(LIBRARIES))
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

VERSIONPY := $(wildcard version.py)

#-- Testsuite
#############
TESTS := first/test.sh second/test.sh

#-- Prefix Testsuite directory name
TESTDIR := Testsuite/
TESTS := $(addprefix $(TESTDIR),$(TESTS))


########################################################################
# TARGETS
########################################################################
# Utility targets (ignore corresponding file names)
.PHONY: all clean $(SUBDIRS) $(SUBCLEAN) prepare_run check run runmpi $(TESTS)

all: $(MODULES)
	$(MAKE) $(SUBDIRS)
	$(MAKE) $(PROGRAMS)
	@echo "MAKE SUCCESSFUL: $(shell date)"

clean: $(SUBCLEAN)
	rm -f *.o *.a *.mod version.inc Makefile.depend $(PROGRAMS)
$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean


prepare_run:
	mkdir $(RUNDIR) 2>/dev/null || rm -f $(RUNDIR)/*
	cd $(RUNDIR) && ln -s $(CURDIR)/Data/* .
	cd $(RUNDIR) && ln -s $(CURDIR)/Input/* .
	cd $(RUNDIR) && ln -s $(CURDIR)/supernu .

run: RUNDIR := $(CURDIR)/Run
run: all prepare_run
	cd $(RUNDIR) && ./supernu

runmpi: RUNDIR := $(CURDIR)/Run
runmpi: all prepare_run
	cd $(RUNDIR) && mpirun -n 2 ./supernu

check: all $(TESTS)
	@echo "TESTSUITE SUCCESSFUL: $(shell date)"
	@echo

########################################################################
# EXPLICIT RULES
########################################################################
#
#-- automatic Makefile generation rule
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
Makefile.compiler:
	$(error USAGE:  cp System/Makefile.compiler.intel-x86_64 Makefile.compiler)

#
#-- program
supernu: supernu.o $(MODULES) $(FILES) $(LIBRARIES)

#
#-- libraries
$(SUBDIRS):
	$(MAKE) -C $@

#
#-- modules
bbxsmod.o bbxsmod.mod: bbxsmod.f elemdatamod.o miscmod.o 
gasgridmod.o gasgridmod.mod: gasgridmod.f90 inputparmod.o 
inputparmod.o inputparmod.mod: inputparmod.f miscmod.o physconstmod.o 
inputstrmod.o inputstrmod.mod: inputstrmod.f elemdatamod.o gasgridmod.o inputparmod.o miscmod.o physconstmod.o 
ionsmod.o ionsmod.mod: ionsmod.f miscmod.o physconstmod.o 
manufacmod.o manufacmod.mod: manufacmod.f gasgridmod.o inputparmod.o miscmod.o physconstmod.o 
miscmod.o miscmod.mod: miscmod.f MISC/lcase.f MISC/warn.f 
mpimod.o mpimod.mod: mpimod.f gasgridmod.o inputparmod.o particlemod.o timestepmod.o timingmod.o
timestepmod.o timestepmod.mod: timestepmod.f90 physconstmod.o 
gammaprofmod.o gammaprofmod.mod: gammaprofmod.f


#
#-- inc files
ifeq ($(VERSIONPY), version.py)
  version.inc: $(VERSIONPY)
	@python2 -B version.inc.py
else
  version.inc: version_dummy.inc
	@cp -vu version_dummy.inc version.inc
endif

#
#-- testsuite
$(TESTS): export RUNDIR = $(CURDIR)/$(patsubst %/,%,$(dir $@))_Run/
$(TESTS):
	$(MAKE) prepare_run
	$(SHELL) $@ $(RUNDIR)
	@echo

#
#-- object dependencies
include Makefile.depend
