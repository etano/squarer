include ../C_OPTIONS

OBJS = cvmgp.o \
       isamax.o isamin.o ismax.o ismin.o \
       mxma.o mxv.o mxva.o \
       sasum.o saxpy.o scopy.o sdot.o ssvdc.o \
       sgeco.o sgedi.o sgefa.o spofa.o sposl.o sscal.o ssum.o sswap.o\
       srotg.o srot.o snrm2.o

forsub : $(OBJS)

clean:
	rm -f *.o core

$(OBJS): mach.p ../C_OPTIONS
