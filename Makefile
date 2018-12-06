#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
include Makefile.include

PROGRAMS := supernu

MODULES := elemdatamod.o nucdatamod.o physconstmod.o randommod.o mpimod.o \
  miscmod.o ionsmod.o ffxsmod.o bfxsmod.o bbxsmod.o \
  inputparmod.o timestepmod.o groupmod.o gridmod.o gasmod.o inputstrmod.o \
  particlemod.o timingmod.o countersmod.o manufacmod.o fluxmod.o totalsmod.o \
  transportmod.o sourcemod.o tbxsmod.o tbsrcmod.o

FILES := banner.o \
  particle_advance.o particle_advance_gamgrey.o fluxtally.o \
  dealloc_all.o read_bbxs_data.o tau_update.o

LIBRARIES := TRANSPORT1/transport1.a TRANSPORT2/transport2.a TRANSPORT3/transport3.a \
  SOURCE/source.a GAS/gas.a GRID/grid.a MISC/misc.a OUTPUT/output.a
SUBDIRS := $(dir $(LIBRARIES))
SUBCLEAN = $(addsuffix .clean, $(SUBDIRS))

SUPERNUDEP := supernu.o $(MODULES) $(FILES) $(LIBRARIES)

VERSIONPY := $(wildcard version.py)
VERSIONDEP := $(filter-out banner.o, $(SUPERNUDEP))
DATE := $(shell date)

#-- Testsuite
#############
TESTS := first/test.sh #second/test.sh

#-- Prefix Testsuite directory name
TESTDIR := Testsuite/
TESTS := $(addprefix $(TESTDIR),$(TESTS))


########################################################################
# TARGETS
########################################################################
# Utility targets (ignore corresponding file names)
.PHONY: all clean $(SUBDIRS) $(SUBCLEAN) prepare_run check run runmpi ready_run $(TESTS) help

all: $(MODULES)
	$(MAKE) $(SUBDIRS)
	$(MAKE) $(PROGRAMS)
	git -C $(dir $(realpath supernu.f90)) diff >gitdiff.txt
	@echo "MAKE SUCCESSFUL: $(shell date)"

clean: $(SUBCLEAN)
	rm -f *.o *.a *.mod version.inc Makefile.dependmod Makefile.depend $(PROGRAMS)
$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean


prepare_run:
	mkdir $(RUNDIR) 2>/dev/null || rm -f $(RUNDIR)/*
	cd $(RUNDIR) && ln -s $(CURDIR)/Data/* .
	cd $(RUNDIR) && ln -s $(CURDIR)/Input/* .
	cd $(RUNDIR) && ln -s $(CURDIR)/supernu .

ready_run: RUNDIR := $(CURDIR)/Run
ready_run: all prepare_run

run: RUNDIR := $(CURDIR)/Run
run: all prepare_run
	cd $(RUNDIR) && ./supernu

runmpi: RUNDIR := $(CURDIR)/Run
runmpi: all prepare_run
	cd $(RUNDIR) && mpirun -n 2 ./supernu

check: all $(TESTS)
	@echo "TESTSUITE SUCCESSFUL: $(shell date)"
	@echo

help:
	@echo
	@echo "Make options: make [help] [all] [clean] [check] [run] [runmpi] [ready_run]"
	@echo
	@echo "help: print this help information"
	@echo "all: compile SuperNu executable"
	@echo "clean: remove executable, object files, and libraries"
	@echo "check: run tests in the Testsuite directory"
	@echo "run: create Run directory and run SuperNu executable in serial"
	@echo "runmpi: create Run directory and run SuperNu executable with MPI on 2 ranks"
	@echo "ready_run: create and set up Run directory, but do not run SuperNu executable"
	@echo

########################################################################
# EXPLICIT RULES
########################################################################
#
#-- automatic Makefile generation rule
Makefile.dependmod:
	$(TOP)/dependmodule.sh $(MODULES) >Makefile.dependmod
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
Makefile.compiler:
	$(error USAGE:  cp System/Makefile.compiler.intel-x86_64 Makefile.compiler)

#
#-- program
supernu: $(SUPERNUDEP)

#
#-- libraries
$(SUBDIRS):
	$(MAKE) -C $@

#
#-- inc files
ifeq ($(VERSIONPY), version.py)
  version.inc: version.py $(VERSIONDEP)
	@python2 -B version.inc.py
	@echo "      build_date = '$(DATE)'" >>$@
else
  version.inc: version_dummy.inc $(VERSIONDEP)
	@cp -v version_dummy.inc version.inc
	@echo "      build_date = '$(DATE)'" >>$@
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
include Makefile.dependmod
include Makefile.depend
