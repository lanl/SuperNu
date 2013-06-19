include Makefile.compiler
include Makefile.inc

########################################################################
# MACROS
########################################################################
PROGRAMS := supernu

MODULES := \
  elemdatamod.o physconstmod.o mpimod.o \
  miscmod.o \
  ionsmod.o ffxsmod.o bfxsmod.o bbxsmod.o \
  inputparmod.o timestepmod.o \
  gasgridmod.o inputstrmod.o \
  particlemod.o \
  timingmod.o

OBJFILES := sourcenumbers.o vacancies.o boundary_source.o interior_source.o particle_advance.o \
 write_output.o diffusion1.o transport1.o read_wlgrid.o wlgrid_setup.o \
 read_bbxs_data.o restart_file.o dealloc_all.o specint.o initialnumbers.o

LIBRARIES := GASGRID/gasgrid.a MISC/misc.a
SUBDIRS := $(dir $(LIBRARIES))

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
bbxsmod.o: physconstmod.o elemdatamod.o miscmod.o
inputparmod.o: miscmod.o physconstmod.o
ionsmod.o: miscmod.o physconstmod.o
miscmod.o: MISC/warn.f MISC/lcase.f
mpimod.o: gasgridmod.o inputparmod.o timestepmod.o timingmod.o particlemod.o
timestepmod.o: physconstmod.o
gasgridmod.o: inputparmod.o
inputstrmod.o: physconstmod.o miscmod.o
#
#-- OBJ FILES
#-- note: prerequisites don't need to include modules as these are always built first
initialnumbers.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o
particle_advance.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o timingmod.o
boundary_source.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o
dealloc_all.o: bbxsmod.o gasgridmod.o ionsmod.o mpimod.o particlemod.o
diffusion1alt.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o
diffusion1.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o
interior_source.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o
read_bbxs_data.o: bbxsmod.o ionsmod.o miscmod.o physconstmod.o timingmod.o
restart_file.o: gasgridmod.o
wlgrid_setup.o: gasgridmod.o physconstmod.o miscmod.o
sourcenumbers.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o
supernu.o: bfxsmod.o ffxsmod.o gasgridmod.o inputparmod.o ionsmod.o mpimod.o particlemod.o physconstmod.o timestepmod.o timingmod.o
transport1.o: gasgridmod.o inputparmod.o particlemod.o physconstmod.o timestepmod.o
vacancies.o: particlemod.o
write_output.o: gasgridmod.o particlemod.o timestepmod.o

banner.o: version.inc mpimod.o
supernu.o: $(MODULES)

#
#-- LIBRARIES
$(SUBDIRS):
	$(MAKE) -C $@
#
#-- TESTSUITE
$(TESTS): export RUNDIR = $(CURDIR)/$(patsubst %/,%,$(dir $@))_Run/
$(TESTS):
	$(MAKE) prepare_run
	$(SHELL) $@ $(RUNDIR)
	@echo
#
#-- PROGRAMS
supernu: $(MODULES) $(OBJFILES) supernu.o banner.o $(LIBRARIES)
