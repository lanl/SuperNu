TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.include

LIBRARY = gasgrid.a
FILES = eos_update.o \
  gasgrid_setup.o gasgrid_update.o gridvolume.o \
  temperature_update.o analytic_source.o \
  analytic_opacity.o physical_opacity.o physical_opacity_subgrid.o \
  leakage_opacity.o fleck_factor.o \
  emission_probability.o analytic_initial.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
