include ../Makefile.compiler
include ../Makefile.inc

########################################################################
# MACROS
########################################################################
LIBRARY = multigroup.a

FILES = analytic_opacity.o analytic_source.o opacity_calculation.o \
  xsections.o

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
analytic_source.o: gasgridmod.mod physconstmod.mod timestepmod.mod particlemod.mod inputparmod.mod
opacity_calculation.o: bbxsmod.mod bfxsmod.mod ffxsmod.mod gasgridmod.mod inputparmod.mod ionsmod.mod miscmod.mod physconstmod.mod timingmod.mod
xsections.o: gasgridmod.mod physconstmod.mod timestepmod.mod