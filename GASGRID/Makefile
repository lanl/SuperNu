include ../Makefile.compiler
include ../Makefile.inc

########################################################################
# MACROS
########################################################################
LIBRARY = gasgrid.a

FILES = eos_update.o opacity_calculation.o \
  gasgrid_setup.o gasgrid_update.o \
  grids.o initials.o xsections.o material_update.o

FILES = grids.o initials.o xsections.o material_update.o

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
gasgrid_setup.o: elemdatamod.mod gasgridmod.mod gstructmod.mod inputparmod.mod miscmod.mod physconstmod.mod rgridmod.mod
gasgrid_update.o: gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod mpimod.mod physconstmod.mod rgridmod.mod timestepmod.mod timingmod.mod
opacity_calculation.o: bbxsmod.mod bfxsmod.mod ffxsmod.mod fluxmod.mod gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod physconstmod.mod timingmod.mod
grids.o: gasgridmod.mod timestepmod.mod
initials.o: physconstmod.mod particlemod.mod gasgridmod.mod timestepmod.mod inputparmod.mod
xsections.o: physconstmod.mod gasgridmod.mod timestepmod.mod
material_update.o: physconstmod.mod gasgridmod.mod timestepmod.mod
