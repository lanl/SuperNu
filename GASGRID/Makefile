include ../Makefile.compiler
include ../Makefile.inc

########################################################################
# MACROS
########################################################################
LIBRARY = gasgrid.a

FILES = eos_update.o opacity_calculation.o \
  gasgrid_setup.o gasgrid_update.o \
  xsections.o material_update.o

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
eos_update.o: ionsmod.mod gasgridmod.mod timestepmod.mod timingmod.mod
gasgrid_setup.o: elemdatamod.mod gasgridmod.mod timestepmod.mod inputparmod.mod miscmod.mod physconstmod.mod
gasgrid_update.o: gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod mpimod.mod physconstmod.mod timestepmod.mod timingmod.mod
opacity_calculation.o: bbxsmod.mod bfxsmod.mod ffxsmod.mod gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod physconstmod.mod timingmod.mod
xsections.o: physconstmod.mod gasgridmod.mod timestepmod.mod
material_update.o: physconstmod.mod gasgridmod.mod timestepmod.mod
