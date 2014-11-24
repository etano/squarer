MACHP   = rtl/mach.p
ENV	= environment

OBJS=asdf

# The compilation options are now in these file
# In all the subdirectories there are includes to these files
#
include C_OPTIONS
include L_OPTIONS
#

default : Nlinpack Nsprng $(ENV) Nrtl Nforsub Nfnlib Nsqdir C_OPTIONS

# machconfig change C_OPTIONS files so that it enforces
#            to recompile all files. 
#            machconfig is executed only if "environment" does not exist
#            The link between rtl and {machine}rtl is now in machconfig
$(ENV) : 
	(machconfig) 

Nlinpack :
	( cd linpack ; sh compile.sh )

Nsprng :
	( cd sprng/SRC ; make )

Nrtl : rtl C_OPTIONS
	( cd rtl ; make rtl )

Nforsub : forsub C_OPTIONS
	( cd forsub ; make forsub )

Nfnlib : fnlib C_OPTIONS
	( cd fnlib ; make fnlib )

Nsqdir : sqdir C_OPTIONS
	( cd sqdir ; make squarer ; make potgen_lr ; make potgen_sr )

clean : 
	(cd linpack; rm *.o *.a)
	(cd sprng/SRC; make clean)
	(cd fnlib; make clean)
	(cd sqdir; make clean)
	(cd forsub; make clean)
	(cd rtl; make clean)
