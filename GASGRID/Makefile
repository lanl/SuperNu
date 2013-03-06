include ../Makefile.compiler
include ../Makefile.inc

########################################################################
# MACROS
########################################################################
LIBRARY = ggrid.a

FILES = eos_update.o opacity_calculation.o \
  ggrid_setup.o ggrid_update.o \

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
eos_update.o: ionsmod.mod ggridmod.mod timestepmod.mod timingmod.mod
ggrid_setup.o: elemdatamod.mod ggridmod.mod gstructmod.mod inputparmod.mod miscmod.mod physconstmod.mod rgridmod.mod
ggrid_update.o: ggridmod.mod inputparmod.mod ionsmod.mod miscmod.mod mpimod.mod physconstmod.mod rgridmod.mod timestepmod.mod timingmod.mod
opacity_calculation.o: bbxsmod.mod bfxsmod.mod ffxsmod.mod fluxmod.mod ggridmod.mod inputparmod.mod ionsmod.mod miscmod.mod physconstmod.mod timingmod.mod
