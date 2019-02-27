#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.include

LIBRARY = grid.a
FILES = grid_setup.o grid_update.o grid_volume.o read_temp_preset.o \
  emitgroup.o emission_probability.o leakage_opacity.o \
  leakage_opacity11.o leakage_opacity1.o leakage_opacity2.o leakage_opacity3.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
