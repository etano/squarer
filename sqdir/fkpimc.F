      program FKPIMC ! revision of 4.29.04
      implicit real*8 (a-h,o-z)
      parameter (mx=127,mdim=3,mdiv=8,mnts=1000,mlevels=8
     +, mterms=10,mnprm=20,mnp=4,mnsys1=1+2**mdiv,md=3)
c if you change this parameter statement, also change those in other subroutines.
#define SIMPLE_SPRNG    
#include "sprng_f.h"
!#define ranf sprng  
      SPRNG_POINTER junkPtr
      integer seed,junk

      dimension r(mx,mdim,mnp,mnsys1),dm(mx,mdiv,4),wt(mx)
     +,intv(2,mnsys1),time(mnsys1),ebeta(mx,mnsys1),spring(mx)
     +,rv(mx),rav(mdim,mnp),virial(mx,mnsys1),cvirial(mx),idim(5)
     +,rdiff(mdim,mnp),u(mdiv,4),uerr(mdiv,4),uextra(4)
     +,exrerr(4),utab(mx,4),av(mdiv,4),sq(mdiv,4),t(mdiv,4),iseed(3)
     +,uexm(4),uexs(4)
c if you change this common statement, also change those in other subroutines.
       common/ch/hbs2m,pot(mnts,4)
     +,dmtable(mnts,4,mterms,mlevels,3),ntermact
       character fname*14,p(20)*28
c************************************************************
c this program calculates the two particle density matrix by
c bisection Monte Carlo
c and compares to that generated in a table.
c for several particles
c************************************************************
      write (*,*)' Begin Feynman-Kacs density matrix computation'
      write (*,*)' Author: D. Ceperley, University of Illinois. Version
     + of 4.29.04'
      open(1,file='fkpimc.in',status='unknown')
      read (1,99) fname
      write (*,*)' Potential and density matrix taken from file ',fname
99    format(a14)
      ln=index(fname,' ')-1
      open (21,file=fname(1:ln)//'.dm',status='unknown')
      write (*,*)'file opened'

      hbs2m=0.d0
      ntypes=0
      ndim=0
1     j=ipickoff(21,p,n,20)
      if(j.eq.1) then
        write (*,*)' data file empty '
        stop
      endif

      if(p(1).eq.'TYPE')then
c this is not correct if particles have different masses
       hbs2m=rlread(p(3))
       ntypes=ntypes+1

      elseif(p(1).eq.'SQUARER') then
        ndim=intread(p(4))
        if(ndim.gt.mdim) stop
        write (*,*)' spatial dimension=',ndim

      elseif(p(1).eq.'GRID') then
        nx=intread(p(2))
        write (*,*) 'GRID ',nx,mnts
        if(nx.le.0.or.nx.gt.mnts)stop
        call setgrid(mnts,nx,rv,p(3))
        write (*,*) 'GRID ',nx,rv(1),rv(nx)
        if(nx.le.0.or.nx.gt.mnts)stop

      endif

      if(ntypes.ne.2.or.ndim.eq.0)go to 1

        write (*,*) ' begin read potential'
        idim(1)=nx
        idim(2)=1
        call checkrnk(21,2,idim)
        call findbegn(21)
        write (*,*)' found begining',nx
        read (21,*) (pot(i,1),i=1,nx)
!       write (6,*)' finished reading potential'
        pot(nx,1)=0.
        call spline(pot,nx,mnts,1)
        write (6,*)' finished reading potential'

c beta is the 1/temperature
      read (1,*) temp
      write (*,*)' temperature=',temp
      beta=1.d0/temp
      read (1,*) ntermact,nd
      write (*,*)' no terms in action ',ntermact
      write (*,*)' number of derivatives ',nd
        do i=1,nd
        call readseln(dmtable(1,1,1,1,i),mnts,4,mterms,mlevels
     +,ntermact,nlvld,beta,nx,i,nx,0)
        enddo

500   close(21)

c point to calculate density matrix
c nparts= number of particles
c nsteps=no. of mc samples to be done (*mx)
c ndiv = number of bisections (work goes as nsteps*2**ndiv)
      read (1,*) nparts,nsteps,ndiv,np
      write  (*,*) ' number particles=',nparts,nsteps,ndiv,np
      nbs=1+2**ndiv
      if(nbs.gt.mnsys1)stop
      if(ndiv.gt.mdiv)stop
      if(nparts.gt.mnp)stop
      if(np.lt.1.or.np.gt.mx)stop

      do ip=1,nparts
         do j=1,2
         read (1,*) (r(1,l,ip,j),l=1,ndim)
         enddo
         write (*,205) ' particle',ip,' initial position'
     +   ,(r(1,l,ip,1),l=1,ndim)
         write (*,205) ' particle',ip,'   final position'
     +   ,(r(1,l,ip,2),l=1,ndim)
205   format(a9,i5,a17,3f12.5)
      enddo
      close(1)
c END OF THE INPUT SECTION

c this is the real-time microsecond clock
!     iseed(1)=irtc()
!     iseed(2)=0
!     iseed(3)=0
c set default random number qseed
!     call iranset(iseed)
      seed = 985456376
      junkPtr = init_sprng(seed,DEFAULT)
      print *, 'Printing information about new stream'
      junk = print_sprng()

c replicate the parallel points
      do 222 k=2,mx
      do 222 ip=1,nparts
      do 222 l=1,ndim
      do 222 id=1,2
222   r(k,l,ip,id)=r(1,l,ip,id)

c initialize final averages
      do 38 k=1,mx
      do 38 l=1,mdiv
      do 38 j=1,4
38    dm(k,l,j)=0.d0

c these are for the virial estimator
      s2=0.d0
      do ip=1,nparts
      do l=1,ndim
      rav(l,ip)=.5d0*(r(1,l,ip,1)+r(1,l,ip,2))
      rdiff(l,ip)=r(1,l,ip,2)-r(1,l,ip,1)
      s2=s2+rdiff(l,ip)**2
      enddo
c     write (*,205) 'rav   part',ip,(rav(l,ip),l=1,ndim)
c     write (*,205) 'rdiff part',ip,(rdiff(l,ip),l=1,ndim)
      enddo
      write (*,*)' diffusion length=',sqrt(2*ndim*hbs2m*beta)
     +,' |R1-R2|=',sqrt(s2)

      time(1)=0.d0
      time(2)=beta
      do 12 j=1,2
12    call en(r(1,1,1,j),mx,mdim,np,ndim,nparts,ebeta(1,j),virial(1,j)
     + ,rav)
      vbar=.5d0*(ebeta(1,1)+ebeta(1,2))
      write (*,*) ' V(R1)=',ebeta(1,1),' V(R2)=',ebeta(1,2)
c     write (*,*) ' virial ',virial(1,1),virial(1,2)

c compute tablulated results
      call toffd(r,r(1,1,1,2),mx,mdim,1,ndim,nparts,utab,nd)
      eprim=nparts*ndim/(2*beta)-s2/(4*hbs2m*beta**2)
      if(md.ge.2)utab(1,2)=utab(1,2)+eprim+vbar
      if(md.ge.3)utab(1,3)=utab(1,3)+eprim
      betai=1.d0/beta
      write (*,*)' number of blocks=',np
      write (*,*)' number of levels=',ndiv
      write (*,*)
      write (*,*)' beginning Monte Carlo, nsteps*nblocks=',nsteps*np

      call second(t0)
c loop over number of monte carlo steps
      do 1001 isteps=1,nsteps
      intv(1,1)=1
      intv(2,1)=2
c loop over number of bisections
      nint=1
      nsl=2
      tau=beta
      do 1000 idiv=1,ndiv
      tau=.5d0*tau
c now bisect each of them
      nint0=nint
      do 100 l=1,nint0
      call bsect(intv,time,l,nint,nsl,hbs2m,r,mx,mdim,mnp,np,ndim
     +   ,nparts)
100   call en(r(1,1,1,nsl),mx,mdim,np,ndim,nparts,ebeta(1,nsl)
     +   ,virial(1,nsl),rav)

      do 960 k=1,np
      spring(k)=.5d0*ndim*nparts/tau
      cvirial(k)=.5d0*ndim*nparts*betai
960   wt(k)=0.d0

      do k=1,nint
      if(intv(1,k).eq.1) k1=intv(2,k)
      if(intv(2,k).eq.2) km=intv(1,k)
      enddo

      cns=0.25d0/(hbs2m*tau*beta)
      do j=1,nparts
      do l=1,ndim
      do k=1,np
      cvirial(k)=cvirial(k)-.5d0*cns*rdiff(l,j)
     +     *(rdiff(l,j)-r(k,l,j,km)+r(k,l,j,k1))
      enddo
      enddo
      enddo

      aa=.5d0/nint
      do l=3,nsl
      do k=1,np
      cvirial(k)=cvirial(k)-aa*virial(k,l)
      enddo
      enddo

      do 965 j=1,nint
      i1=intv(1,j)
      i2=intv(2,j)
      do 965 i=1,nparts
      do 965 l=1,ndim
      do 965 k=1,np
965   spring(k)=spring(k)-cns*(r(k,l,i,i1)-r(k,l,i,i2))**2

c integrate them up
      do 975 l=1,nsl
      ww=tau
      if(l.le.2) ww=.5d0*tau
      do 975 k=1,np
975   wt(k)=wt(k)+ww*ebeta(k,l)

      do k=1,np
      w=exp(-wt(k))
      dm(k,idiv,1)=dm(k,idiv,1)+w
      dm(k,idiv,2)=dm(k,idiv,2)+w*(wt(k)*betai+cvirial(k))
      dm(k,idiv,3)=dm(k,idiv,3)+w*cvirial(k)
      dm(k,idiv,4)=dm(k,idiv,4)+w*(spring(k)-cvirial(k))
      enddo


1000  continue
c write out a path for visualization
c     if(ifout.eq.0) then
c     do k=1,np
c     if(ranf().lt.exp(vcut-wt(k))) then
c        ifout=1
c write out configuration in normal time order
c         tnow=0
c         tau=beta/(nsl-1)
c         do i1=1,nsl
c find minimum t
c            do i2= 1,nsl
c            if(abs(time(i2)-tnow).le..1*tau) then
c              write (20,21) time(i2),(r(k,l,i2),l=1,ndim)
21             format(4f10.5)
c            endif
c            enddo
c         tnow=tnow+tau
c         enddo
c         go to 1001
c     endif
c     enddo
c     endif
 1001  continue

      write (*,*)'level  Action      error       E virial    error     
     + T virial    error'
      write (*,1020)0,beta*vbar,0.,eprim+vbar,0.,eprim,0.
c find the mean and standard deviation of dm and energy
      nblock=np
c ASSUMING tau^2 extrapolation
      n=2
      wex=1.d0/(2.0d0**n-1.d0)
      wex1=wex+1.d0

      do l=1,4
          uexm(l)=0.d0
          uexs(l)=0.d0
          do idiv=1,ndiv
            av(idiv,l)=0.d0
            sq(idiv,l)=0.d0
            t(idiv,l)=0.d0
          enddo
      enddo

      do k=1,nblock
          do idiv=1,ndiv
          dm1=dm(k,idiv,1)
          do l=1,4
c take log or normalize
           if(l.eq.1) then
              value=-log(dm1/(nsteps))
           else
              value=dm(k,idiv,l)/dm1
           endif
           av(idiv,l) =av(idiv,l) +value
           sq(idiv,l) =sq(idiv,l) +value**2
           t(idiv,l)=t(idiv,l)+dm(k,idiv,l)
           dm(k,idiv,l)=value
         enddo
         enddo
c extrapolate
         do l=1,4
           uex=wex1*dm(k,ndiv,l)-wex*dm(k,ndiv-1,l)
           uexm(l)=uexm(l)+uex
           uexs(l)=uexs(l)+uex**2
         enddo
      enddo

      do l=1,4
        do idiv=1,ndiv
         av(idiv,l)=av(idiv,l)/nblock
         uerr(idiv,l)=sqrt((sq(idiv,l)/nblock-av(idiv,l)**2)
     +      /(nblock-1.d0))
         if(l.eq.1) then
            u(idiv,1)=-log(t(idiv,1)/(nsteps*nblock))
         else
            u(idiv,l)=t(idiv,l)/t(idiv,1)
         endif
        enddo
         uextra(l)=wex1*u(ndiv,l)-wex*u(ndiv-1,l)
         exrerr(l)=sqrt((uexs(l)/nblock-(uexm(l)/nblock)**2)
     +     /(nblock-1.d0))
      enddo

      do idiv=1,ndiv
      write (*,1020) idiv,(u(idiv,l),uerr(idiv,l),l=1,nd)
      enddo
      write (*,1022) (uextra(l),exrerr(l),l=1,nd)
      write (*,1024) (utab(1,l),l=1,nd)
      write (*,1021) 'MC-tab',((uextra(l)-utab(1,l)),l=1,nd)
      write (*,1021) 'chi   ',((uextra(l)-utab(1,l))/exrerr(l),l=1,nd)
      write (*,1021) 'bias  ',((u(ndiv,l)-av(ndiv,l))/nblock,l=1,nd)
      write (*,1021) 'rel.ex',((uextra(l)-u(ndiv,l))/exrerr(l),l=1,nd)
1020  format('     ',i1,6e12.4)
1021  format(a6,3(12x,e12.4))
1022  format('extrp ',8e12.4)
1024  format('table ',4(e12.4,12x))

      if(nd.ge.2) then
      do idiv=1,ndiv
      write (*,1025) idiv, u(idiv,4),uerr(idiv,4)
     +,u(idiv,4)/uerr(idiv,4)
1025  format(' virial test:level=',i3,' e_T-e_V=',e12.4,' error=',e12.4
     +,' ratio=',f9.4)
      enddo
      endif

      call second(t1)
      write (*,1040) t1-t0
1040  format(/' time ',f20.3,'CPU seconds')
      end
 
      subroutine en(r,mx,mdim,nx,ndim,nparts,ebeta,f,rbar)
      implicit real*8 (a-h,o-z)
      parameter (mxx=127,mnts=1000,mlevels=8,mterms=10)
       common/ch/hbs2m,pot(mnts,4)
     +,dmtable(mnts,4,mterms,mlevels,3),ntermact
      dimension r(mx,mdim,nparts),ebeta(mx),rr(mxx),p(mxx),ind(mxx)
     +,f(mx),dvdr(mxx),rbar(mdim,nparts),rdpdr(mxx)

      do k=1,nx
      f(k)=0.d0
      ebeta(k)=0.d0
      enddo

      do i=2,nparts
      do j=1,i-1

      if(ndim.eq.1) then
      do 1 k=1,nx
1     rr(k)=r(k,1,i)-r(k,1,j)

      elseif(ndim.eq.2) then
      do 2 k=1,nx
2     rr(k)=sqrt((r(k,1,i)-r(k,1,j))**2 +(r(k,2,i)-r(k,2,j))**2)

      elseif(ndim.eq.3) then
      do 3 k=1,nx
3     rr(k)=sqrt( (r(k,1,i)-r(k,1,j))**2 +(r(k,2,i)-r(k,2,j))**2
     + + (r(k,3,i)-r(k,3,j))**2   )
      endif

c convert to a grid
      call ti(rr,ind,p,nx,rdpdr)
      call evpol(pot,ind,p,nx,mnts,ebeta,rdpdr,dvdr)

      do l=1,ndim
      do k=1,nx
      f(k)=f(k)-dvdr(k)* (r(k,l,i)-r(k,l,j) )
     +   *(r(k,l,i)-r(k,l,j)-rbar(l,i)+rbar(l,j))
      enddo
      enddo

      enddo
      enddo

      return
      end

      subroutine ti(r,ind,rem,n,rdxdr)
      implicit real*8 (a-h,o-z)
      dimension ind(n),r(n),rem(n),rdxdr(n)
      character gtype*8
      common/gcom/gtype,dri,di,x0,rr,nx

      if(gtype.eq.'LINEAR') then
      do j=1,n
      x=x0+dri*r(j)
      ind(j)=max(1,min(nx,int(x)))
      rem(j)=x-ind(j)
      rdxdr(j)=dri/r(j)
      enddo

      elseif(gtype.eq.'POWER') then
      do j=1,n
      x=x0+(dri*(r(j)-rr))**di
      ind(j)=max(1,min(nx,int(x)))
      rem(j)=x-ind(j)
      rdxdr(j)=0.
      enddo
      stop

      elseif(gtype.eq.'LOG') then
      do j=1,n
      x=1+x0*log(dri*r(j))
      ind(j)=max(1,min(nx,int(x)))
      rem(j)=x-ind(j)
      rdxdr(j)=x0/r(j)**2
      enddo

      elseif(gtype.eq.'INVERSE') then

      do j=1,n
          x=x0+dri/(r(j)-rr)
          ind(j)=max(1,min(nx,int(x)))
          rem(j)=x-ind(j)
          rdxdr(j)=- dri/((r(j)-rr)**2*r(j))
      enddo

      endif

      return
      end

      subroutine evpol(pot,ind,p,mx,mnts,ebeta,rdpdr,dvdr)
      implicit real*8 (a-h,o-z)
      dimension pot(mnts,4),ebeta(mx),ind(mx),p(mx),dvdr(mx),rdpdr(mx)
      do 4 k=1,mx
      dvdr(k)=rdpdr(k)*(pot(ind(k),2)+p(k)*(2.d0*pot(ind(k),3)
     +    +3.d0*p(k)*pot(ind(k),4) ) )
4     ebeta(k)=ebeta(k)+pot(ind(k),1)+p(k)*(pot(ind(k),2)+p(k)*
     +  (pot(ind(k),3)+p(k)*pot(ind(k),4))   )
      return
      end

      subroutine toffd(r1,r2,mx,mdim,nx,ndim,nparts,u,nd)
      implicit real*8 (a-h,o-z)
      parameter (mxx=127,mnts=1000,mlevels=8,mterms=10)
       common/ch/hbs2m,pot(mnts,4)
     +,dmtable(mnts,4,mterms,mlevels,3),ntermact
      dimension r1(mx,mdim,nparts),r2(mx,mdim,nparts),u(mx,3)
     + ,rr1(mxx),rr2(mxx),p(mxx),ind(mxx),rb(mxx)
     +,dvdr(mxx),rdpdr(mxx), var(mxx,10)
     +,uend(mxx,3)
 
c we will calculate u and udot
      do l=1,nd
      do k=1,nx
      uend(k,l)=0.d0
      u(k,l)=0.d0
      enddo 
      enddo 
 
      do i=2,nparts
      do j=1,i-1
 
      if(ndim.eq.1) then
      do 1 k=1,nx
      var(k,1)=(r1(k,1,i)-r1(k,1,j)-r2(k,1,i)+r2(k,1,j))**2 
      rr1(k)=r1(k,1,i)-r1(k,1,j)
1     rr2(k)=r2(k,1,i)-r2(k,1,j)
 
      elseif(ndim.eq.2) then
      do 2 k=1,nx
      var(k,1)=(r1(k,1,i)-r1(k,1,j)-r2(k,1,i)+r2(k,1,j))**2 
     +      +  (r1(k,2,i)-r1(k,2,j)-r2(k,2,i)+r2(k,2,j))**2 
      rr1(k)=sqrt((r1(k,1,i)-r1(k,1,j))**2 +(r1(k,2,i)-r1(k,2,j))**2)
2     rr2(k)=sqrt((r2(k,1,i)-r2(k,1,j))**2 +(r2(k,2,i)-r2(k,2,j))**2)
 
      elseif(ndim.eq.3) then
      do 3 k=1,nx
      var(k,1)=(r1(k,1,i)-r1(k,1,j)-r2(k,1,i)+r2(k,1,j))**2 
     +    +    (r1(k,2,i)-r1(k,2,j)-r2(k,2,i)+r2(k,2,j))**2 
     +    +    (r1(k,3,i)-r1(k,3,j)-r2(k,3,i)+r2(k,3,j))**2 
      rr1(k)=sqrt( (r1(k,1,i)-r1(k,1,j))**2 +(r1(k,2,i)-r1(k,2,j))**2
     + + (r1(k,3,i)-r1(k,3,j))**2   )
3     rr2(k)=sqrt( (r2(k,1,i)-r2(k,1,j))**2 +(r2(k,2,i)-r2(k,2,j))**2
     + + (r2(k,3,i)-r2(k,3,j))**2   )
      endif
 
c the end point term for r1
      call ti(rr1,ind,p,nx,rdpdr)
      do l=1,nd
      call evpol(dmtable(1,1,1,1,l),ind,p,nx,mnts,uend(1,l),rdpdr,dvdr)
      enddo
c the end point term for r2
      call ti(rr2,ind,p,nx,rdpdr)

      do l=1,nd
      call evpol(dmtable(1,1,1,1,l),ind,p,nx,mnts,uend(1,l),rdpdr,dvdr)
      enddo
 

c now off-diagonal terms
      if(ntermact.gt.1) then
        do k=1,nx
        rb(k)=.5d0*(rr1(k)+rr2(k))
        var(k,2)=(rr1(k)-rr2(k))**2
        enddo

        call ti(rb,ind,p,nx,rdpdr)

       do l=1,nd
       call vsum2(nx,ind,p,dmtable(1,1,2,1,l),mnts,var,ntermact-1
     +  ,mxx,4,ndim,u(1,l))
        enddo
      endif
 
      enddo
      enddo
        do l=1,nd
        do k=1,nx
        u(k,l)=u(k,l)+.5d0*uend(k,l)
        enddo
        enddo
 
      return
      end

      subroutine bsect(intv,time,l,nint,nsl,hbs2m,r,mx,mdim,mnp,nx,ndim
     +   ,nparts)
      implicit real*8 (a-h,o-z)
      dimension intv(2,1),time(1),r(mx,mdim,mnp,1)
       i1=intv(1,l)
       i2=intv(2,l)
      nint=nint+1
      nsl=nsl+1
      time(nsl)=.5d0*(time(i1)+time(i2))
      width=hbs2m*(time(i2)-time(i1))
      if(width.le.0)then
        write (*,*)' width<0 ',width,time(i2),time(i1),hbs2m
     & , nint,nsl,nx,nparts,i1,i2,(time(l),l=1,nsl)
        write (*,*) l,(intv(1,k),intv(2,k),k=1,l)
        stop
      endif
      intv(2,l)=nsl
      intv(1,nint)=nsl
      intv(2,nint)=i2
      do 100 ip=1,nparts
      do 100 j=1,ndim
      do 100 k=1,nx
100   r(k,j,ip,nsl)=cos(3.14159265357d0*sprng())
     +*sqrt(-width*log(sprng()))
     ++.5d0*(r(k,j,ip,i1)+r(k,j,ip,i2))
      return
      end

      subroutine vsum2(m,i,p,t,n,var,nterm,mnp,morder,idim,u)
      implicit real*8 (a-h,o-z)
c we assume the ordering of terms from squarer
c idim=1 1 d tables , idim=2 full tables.

      real*8 p,t,var
      integer i,n,nterm,mnp,morder,nkt,k,l,j,m,idim

      dimension i(m),p(m),t(n,morder,nterm),var(mnp,nterm),u(m)

c determine the expansion coeeficents

      if(idim.eq.1) then
         do 400 l=2,nterm
         do 400 j=1,m
400      var(j,l)=var(j,l-1)*var(j,1)
      else
       if(nterm.gt.2) then
        nkt=2
        k=1
c go to one higher order
50      k=k+1
c multiply by s**2
        do 30 l=1,k
        do 30 j=1,m
30      var(j,nkt+l)=var(j,nkt+l-k)*var(j,1)
c multiply by z**2
        do j=1,m
          var(j,nkt+k+1)=var(j,nkt)*var(j,2)
        enddo

        nkt=nkt+k+1
        if(nkt.gt.nterm) then
           write (*,*) 'problem in vsum2 ',nterm,nkt
           stop
        elseif(nkt.lt.nterm) then
             go to 50
        endif

        endif
      endif

c now do table look-up-for nterm large the loops should be inverted
      do 40 l=1,nterm
        do 40 j=1,m
40       u(j)=u(j)+(t(i(j),1,l)
     &                 +p(j)*(t(i(j),2,l)
     &                        +p(j)*(t(i(j),3,l)
     &                               +p(j)*t(i(j),4,l))))*var(j,l)
      return
      end

