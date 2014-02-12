ifneq ($(TOP), $(CURDIR))
%.mod: $(TOP)/%.mod
	ln -f $< $(CURDIR)/
endif

%.inc: $(TOP)/%.inc
	ln -f $< $(CURDIR)/

%.o: %.f
	$(FC) $(FCFLAGS) $(FCFLAGS_FF) -c $<
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.a: $(FILES) $(MODULES)
	$(AR) $(ARFLAGS) $@ $^
#	$(RANLIB) $@
