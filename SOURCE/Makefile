TOP := $(CURDIR)/..

include $(TOP)/Makefile.compiler
include $(TOP)/Makefile.include

LIBRARY = source.a
FILES := sourceenergy.o sourceenergy_analytic.o sourceenergy_misc.o \
  sourcenumbers.o vacancies.o \
  boundary_source.o interior_source.o initialnumbers.o initial_particles.o \
  source_transformdirection.o analytic_initial.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
