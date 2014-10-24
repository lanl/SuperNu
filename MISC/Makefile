TOP := $(CURDIR)/..

include ../Makefile.compiler
include ../Makefile.include

LIBRARY = misc.a
FILES = locate.o sorti.o memusg.o
#-- the following routines are included in miscmod.f in order to avoid
#-- explicit interfaces:
# warn.o lcase.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
