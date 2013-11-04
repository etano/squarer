include OPTIONS

OBJS = i1mach.o irtc.o pack.o r1mach.o second.o flush.o	unpack.o timedate.o

default:
	echo "Type 'make' in the main fnpimc directory."

rtl: mach.p $(OBJS)
  
clean:
	rm -f *.o core

$(OBJS) : mach.p OPTIONS
