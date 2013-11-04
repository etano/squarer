include OPTIONS

# iranget.o iranset.o ranf.o
# i1mach.o r1mach.o second.o 
# pack.o unpack.o 
OBJS = i1mach.o irtc.o pack.o r1mach.o second.o timedate.o unpack.o

default:
	echo "Type 'make' in the main fnpimc directory."

rtl: mach.p $(OBJS)
  
clean:
	rm -f *.o core a.out

$(OBJS) : mach.p OPTIONS
