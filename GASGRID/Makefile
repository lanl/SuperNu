include ../Makefile.compiler
include ../Makefile.inc

########################################################################
# MACROS
########################################################################
LIBRARY = gasgrid.a

FILES = analytic_opacity.o eos_update.o opacity_calculation.o \
  gasgrid_setup.o gasgrid_update.o gasgrid_upold.o\
  xsections.o temperature_update.o analytic_source.o

########################################################################
# TARGETS
########################################################################
# Utility targets
.PHONY: clean

$(LIBRARY): $(FILES)
	$(AR) lruv $@ $^
	$(RANLIB) $@

clean:
	rm -f *.o *.mod *.MOD $(LIBRARY)

########################################################################
# EXPLICIT RULES
########################################################################
analytic_opacity.o: gasgridmod.mod physconstmod.mod
analytic_source.o: gasgridmod.mod physconstmod.mod
eos_update.o: gasgridmod.mod ionsmod.mod timestepmod.mod timingmod.mod
gasgrid_setup.o: elemdatamod.mod gasgridmod.mod inputparmod.mod miscmod.mod physconstmod.mod timestepmod.mod
gasgrid_update.o: gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod mpimod.mod physconstmod.mod timestepmod.mod timingmod.mod
##
gasgrid_upold.o: gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod mpimod.mod physconstmod.mod timestepmod.mod timingmod.mod
##
opacity_calculation.o: bbxsmod.mod bfxsmod.mod ffxsmod.mod gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod physconstmod.mod timingmod.mod
temperature_update.o: gasgridmod.mod physconstmod.mod timestepmod.mod
xsections.o: gasgridmod.mod physconstmod.mod timestepmod.mod
