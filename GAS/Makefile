TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.include

LIBRARY = gas.a
FILES = eos_update.o \
  gas_setup.o gas_update.o \
  temperature_update.o analytic_source.o \
  analytic_opacity.o physical_opacity.o physical_opacity_subgrid.o \
  leakage_opacity1.o leakage_opacity2.o leakage_opacity3.o \
  fleck_factor.o opacity_planckmean.o \
  emission_probability.o analytic_initial.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
