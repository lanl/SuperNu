include ../Makefile.compiler
include ../Makefile.inc

########################################################################
# MACROS
########################################################################
LIBRARY = gasgrid.a

FILES = eos_update.o \
  gasgrid_setup.o gasgrid_update.o \
  temperature_update.o analytic_source.o \
  analytic_opacity.o physical_opacity.o convert_cap2capros.o \
  leakage_opacity.o fleck_factor.o \
  emission_probability.o

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
analytic_source.o: gasgridmod.mod physconstmod.mod timestepmod.mod particlemod.mod inputparmod.mod manufacmod.mod
eos_update.o: gasgridmod.mod ionsmod.mod timestepmod.mod timingmod.mod
gasgrid_setup.o: elemdatamod.mod gasgridmod.mod inputparmod.mod miscmod.mod physconstmod.mod timestepmod.mod inputstrmod.mod manufacmod.mod
gasgrid_update.o: gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod mpimod.mod physconstmod.mod timestepmod.mod timingmod.mod
##
physical_opacity.o: bbxsmod.mod bfxsmod.mod ffxsmod.mod gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod physconstmod.mod timingmod.mod
convert_cap2capros.o: gasgridmod.mod
temperature_update.o: gasgridmod.mod physconstmod.mod timestepmod.mod manufacmod.mod
leakage_opacity.o: gasgridmod.mod timestepmod.mod particlemod.mod physconstmod.mod
emission_probability.o: gasgridmod.mod physconstmod.mod
fleck_factor.o: gasgridmod.mod timestepmod.mod physconstmod.mod
