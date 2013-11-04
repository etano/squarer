include ../C_OPTIONS
include ../L_OPTIONS

LIBFN=   ../fnlib/*.o
LIBRTL=  ../rtl/*.o
LIBFOR=  ../forsub/*.o

OBJS =  denmat.o strans.o vmsbi.o fdenmt.o isym.o vcos.o primitive.o \
        initialu.o vlegp.o sqit.o diagonal.o pwc.o \
	smear.o fitter.o sqit_linear.o epp.o breakup.o

squarer : squarer.o $(OBJS)
	$(FC) $@.o $(LFLAGS) $(OBJS) $(LIBFN)   \
	$(LIBFOR) $(LIBRTL) $(LDFLAGS) -o $@

linear_squarer : linear_squarer.o $(OBJS)
	$(FC) $@.o $(LFLAGS) $(OBJS) $(LIBFN)   \
	$(LIBFOR) $(LIBRTL) $(LDFLAGS) -o $@

potgen_sr : potgen_sr.o $(OBJS)
	$(FC) $@.o $(LFLAGS) $(OBJS)  $(LIBFOR)  \
	$(LIBFN) $(LIBRTL) $(LDFLAGS)  -o $@

potgen_lr : potgen_lr.o $(OBJS)
	$(FC) $@.o $(LFLAGS) $(OBJS)  $(LIBFOR)  \
	$(LIBFN) $(LIBRTL) $(LDFLAGS)  -o $@

nlook : nlook.o
	$(FC) $@.o $(LFLAGS) $(LIBDIR)   $(LIBFOR) \
	$(LIBFN) $(LIBRTL) $(LDFLAGS)  -o $@

diffdm : diffdm.o
	$(FC) $@.o $(LFLAGS) $(LIBDIR)   $(LIBFOR) \
	$(LIBRTL)   -o $@

fkpimc : fkpimc.o
	$(FC) $@.o $(LFLAGS) $(LIBDIR)   $(LIBFOR) $(OBJS) $(LIBFN)\
	$(LIBRTL) $(LDFLAGS)  -o $@ 

fkpimc.o : fkpimc.F
	$(FC) $(FFLAGS) $(LIBRTL) $(LDFLAGS) -c fkpimc.F

tidy:
	- rm -f core *.o *.dvi *.ps *.aux *.log

clean: tidy
	- rm -f squarer fkpimc nlook potgen_lr potgen_sr diffdm 

$(OBJS):  ../C_OPTIONS
