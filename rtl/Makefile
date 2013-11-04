include ../C_OPTIONS

# iranget.o iranset.o ranf.o
# i1mach.o r1mach.o second.o 
# pack.o unpack.o 
OBJS =  r1mach.o i1mach.o pack.o timedate.o unpack.o second.o flush.o

default:
	echo "Type 'make' in the main fnpimc directory."

rtl: mach.p $(OBJS)
  
clean:
	rm -f *.o core 

$(OBJS) : mach.p OPTIONS
