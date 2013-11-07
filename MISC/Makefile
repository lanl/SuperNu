TOP := $(CURDIR)/..

include ../Makefile.compiler
include ../Makefile.inc

FILES = locate.o sorti.o memstats.o
#-- the following routines are included in miscmod.f in order to avoid
#-- explicit interfaces:
# warn.o lcase.o

.PHONY: all clean

all: misc.a
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

misc.a: $(FILES)

include Makefile.depend
