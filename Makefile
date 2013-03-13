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
 read_bbxs_data.o restart_file.o dealloc_all.o

LIBRARIES = GASGRID/gasgrid.a MISC/misc.a
SUBDIRS = $(dir $(LIBRARIES))

VERSIONPY = $(wildcard version.py)

#-- Testsuite
########################################################################
TESTS := first/test.sh second/test.sh

#-- Prefix Testsuite directory name
TESTDIR := Testsuite/
TESTS := $(addprefix $(TESTDIR),$(TESTS))


########################################################################
#-- CONSTANTS
########################################################################
RUNDIR = $(CURDIR)/Run

########################################################################
# TARGETS
########################################################################
# Utility targets (ignore corresponding file names)
.PHONY: all put-intrepid clean check run tar $(SUBDIRS) prepare_run $(TESTS)

all: $(MODULES) $(SUBDIRS) $(PROGRAMS)
	@echo "MAKE SUCCESSFUL: $(shell date)"
	@echo

clean:
	for d in $(SUBDIRS); do ($(MAKE) -C $$d clean); done
	rm -f *.o *.a *.mod *.MOD version.inc
	rm -f $(PROGRAMS)

prepare_run:
	mkdir $(RUNDIR) 2>/dev/null || rm -f $(RUNDIR)/*
	cd $(RUNDIR) && ln -s $(CURDIR)/Data/* .
	cd $(RUNDIR) && ln -s $(CURDIR)/Input/* .
	cd $(RUNDIR) && ln -s $(CURDIR)/supernu .

run: all prepare_run
	cd $(RUNDIR) && ./supernu

check: all $(TESTS)
	@echo "TESTSUITE SUCCESSFUL: $(shell date)"
	@echo


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
#-- TESTSUITE
$(TESTS):
	$(MAKE) prepare_run
	$(SHELL) $@
	@echo
#
#-- PROGRAMS
supernu: $(MODULES) $(OBJFILES) supernu.o banner.o $(LIBRARIES)
