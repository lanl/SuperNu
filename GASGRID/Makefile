TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.inc

FILES = eos_update.o \
  gasgrid_setup.o gasgrid_update.o \
  temperature_update.o analytic_source.o \
  analytic_opacity.o physical_opacity.o convert_cap2capros.o \
  leakage_opacity.o fleck_factor.o \
  emission_probability.o analytic_initial.o

.PHONY: all clean

all: gasgrid.a
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD $(LIBRARY)

gasgrid.a: $(FILES)

include Makefile.depend
