include ../Makefile.compiler
include ../Makefile.inc

########################################################################
# MACROS
########################################################################
LIBRARY = misc.a

FILES = locate.o sorti.o memstats.o \
#-- the following routines are included in miscmod.f in order to avoid
#-- explicit interfaces:
# warn.o lcase.o

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
memstats.o: miscmod.mod
