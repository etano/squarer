include ../C_OPTIONS
include ../L_OPTIONS


OBJS = alngam.o asyik.o alpset.o \
       besi.o bessi.o bessi0.o bessj0.o bi1.o \
       checkin.o checkm.o checkrnk.o chgclm.o chgrow.o convert.o \
       cossin.o csevl.o \
       deltar.o dleak.o pairdnode.o dnode.o dtabl1.o \
       distable.o delta_action.o\
       echo.o eupdat.o \
       fillk.o findbegn.o fitpn.o flip.o fpke.o \
       gamma.o gamlim.o gammi.o ggrid.o gradtabl.o \
       ifind.o inits.o interp.o intread.o invert.o ipickoff.o ipwrt.o \
       locate.o lsapr.o \
       mcheck.o mdlng.o midpoint.o mmbsjn.o multiup.o \
       newdet.o newup.o nflip.o ntime.o \
       plint2.o plint3.o pmult.o prspline.o \
       r9lgmc.o readcon.o readseln.o rlread.o rmin.o \
       samcg.o samfr.o sang.o setav.o setfpdm.o setgrid.o \
       shells.o sites.o spline.o sumavv.o symm.o \
       tabl.o upack.o svdcmp.o svbksb.o \
       va10a.o vandvir.o vsum0.o vsum1.o \
       wpack.o writecon.o \
      zeroav.o mnbrak.o brent.o  dbrent.o \
       gvsum.o isrchfgt.o indexx.o fitpnnew.o splint2D.o splint3D.o basis.o \
	computespl2.o struve.o bessj1.o zpolint3.o
fnlib: $(OBJS)

clean:
	rm -f *.o core 

$(OBJS): mach.p ../C_OPTIONS

deltar.o distable.o midpoint.o rmin.o samcg.o samfr.o: cpbc.cm

setav.o sumavv.o zeroav.o: caver.cm

.F.o :
	$(FC) -c $*.F  $(LDFLAGS)
.SUFFIXES : .o .F

