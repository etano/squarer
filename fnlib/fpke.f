      subroutine fpke(nppss,rknorm,kmult,vol,hbs2m,enorm,nspins,ndim
     +               ,energy,mnsh)
      include 'mach.p'
      integer nppss,kmult,nspins,ndim,i,needed,mult,kl,ks,mnsh
      real*8 rknorm,vol,enorm,energy,pi,sangle,tktot,tkitot,tkf 
      real*8 fermiwv,tkinf,hbs2m
      dimension rknorm(mnsh),kmult(0:mnsh),nppss(nspins)
c determine the infinite system and finite system free particle kinetic energies
      pi=-1.
      pi=acos(pi)
c nppss=number of states occupied. rknorm=list of k magnitudes
c kmult are the multiplicites of shells (as order by shell)
c vol is the volume of the box. hbs2m=hbar**2/2*mass, enorm is the
c energy conversion. nspins=number of spin states, ndim= dimensionality
      sangle=2*pi*(ndim-1)
      tktot=0.
      tkitot=0.
      do 1 i=1,nspins
       tkf=0.
c fill up one at the origin
       needed=nppss(i)-1
       do 2 ks=1,3000
        if(ks.eq.1) mult=2*kmult(1)
        if(ks.gt.1) mult=2*(kmult(ks)-kmult(ks-1))
        if(mult.le.0) then
	 write(6,*) 'stop in fpke, mult<0 ',mult
	 stop
        end if
        mult=min0(needed,mult)
        tkf=tkf+hbs2m*rknorm(ks)**2*mult
        needed=needed-mult
        kl=ks
c     write (6,*) ' shell ',ks,' mult ',mult,' k ',rknorm(ks)
        if(needed.le.0) go to 3
2      continue
       write (*,*)' too few states in fpke '
       stop
3      continue
       fermiwv=2*pi*(ndim*nppss(i)/(vol*sangle))**(1./ndim)
       tkinf=hbs2m*sangle*vol*fermiwv**(ndim+2)/((2.*pi)**ndim*(ndim+2))
       write (6,*)'  spin ',i,' ke finite ',tkf,' infinite ',tkinf
       write (6,*) ' fermiwv ',kl,rknorm(kl),fermiwv
       tktot=tktot+tkf
       tkitot=tkitot+tkinf
1     continue
      tktot=tktot/enorm
      tkitot=tkitot/enorm
      write (6,*)' free particle ke finite ',tktot,' infinite ',tkitot
     +,' difference ',tktot-tkitot
      energy=enorm*tktot
      return
      end
