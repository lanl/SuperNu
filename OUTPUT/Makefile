#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.include

LIBRARY = output.a
FILES = open_logfiles.o output_grid.o output_flux.o output_gamflux.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
