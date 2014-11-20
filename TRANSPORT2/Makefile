TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.include

LIBRARY = transport2.a
FILES = diffusion2.o transport2.o advection2.o transport2_gamgrey.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
