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

MODULES := elemdatamod.o physconstmod.o mpimod.o \
  miscmod.o ionsmod.o ffxsmod.o bfxsmod.o bbxsmod.o \
  inputparmod.o timestepmod.o groupmod.o gridmod.o gasmod.o inputstrmod.o \
  particlemod.o timingmod.o manufacmod.o fluxmod.o totalsmod.o

FILES := sourceenergy.o sourceenergy_analytic.o sourceenergy_gamma.o \
  sourcenumbers.o vacancies.o \
  boundary_source.o interior_source.o particle_advance.o \
  particle_advance_gamgrey.o write_output.o wlgrid_setup.o \
  read_bbxs_data.o restart_file.o dealloc_all.o initialnumbers.o \
  initialnumbers.o initial_particles.o tau_update.o \
  banner.o

LIBRARIES := TRANSPORT1/transport1.a TRANSPORT2/transport2.a TRANSPORT3/transport3.a \
  GAS/gas.a  GRID/grid.a MISC/misc.a
SUBDIRS := $(dir $(LIBRARIES))
SUBCLEAN = $(addsuffix .clean, $(SUBDIRS))

SUPERNUDEP := supernu.o $(MODULES) $(FILES) $(LIBRARIES)

VERSIONPY := $(wildcard version.py)
VERSIONDEP := $(filter-out banner.o, $(SUPERNUDEP))
DATE := $(shell date)

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
.PHONY: all clean $(SUBDIRS) $(SUBCLEAN) prepare_run check run runmpi ready_run $(TESTS)

all: $(MODULES)
	$(MAKE) $(SUBDIRS)
	$(MAKE) $(PROGRAMS)
	hg --cwd $(dir $(realpath supernu.f90)) diff >hgdiff.txt
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

##
##-- modules
#bbxsmod.o bbxsmod.mod: bbxsmod.f elemdatamod.o miscmod.o 
#gasmod.o gasmod.mod: gasmod.f90 inputparmod.o 
#inputparmod.o inputparmod.mod: inputparmod.f miscmod.o physconstmod.o 
#inputstrmod.o inputstrmod.mod: inputstrmod.f elemdatamod.o gasmod.o inputparmod.o miscmod.o physconstmod.o 
#ionsmod.o ionsmod.mod: ionsmod.f miscmod.o physconstmod.o 
#manufacmod.o manufacmod.mod: manufacmod.f gasmod.o inputparmod.o miscmod.o physconstmod.o 
#miscmod.o miscmod.mod: miscmod.f MISC/lcase.f MISC/warn.f 
#mpimod.o mpimod.mod: mpimod.f gasmod.o inputparmod.o particlemod.o timestepmod.o timingmod.o
#timestepmod.o timestepmod.mod: timestepmod.f90 physconstmod.o 
#profiledatamod.o profiledatamod.mod: profiledatamod.f


#
#-- inc files
ifeq ($(VERSIONPY), version.py)
  version.inc: version.py $(VERSIONDEP)
	@python2 -B version.inc.py
	@echo "      data build_date /'$(DATE)'/" >>$@
else
  version.inc: version_dummy.inc $(VERSIONDEP)
	@cp -vu version_dummy.inc version.inc
	@echo "      data build_date /'$(DATE)'/" >>$@
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
