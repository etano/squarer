      program potgen
      implicit real*8 (a-h,o-z)
c generates a potential table Updated to include optimized coulomb 4/30/93
      parameter (mn=10,mx=5000)
      dimension pott(mx,2),rv(mx),diff(mx),expon(20),cexpon(20)
      dimension ypot(mx),hbs2m(2),rk(mx),potlrk(mx)
      character fname*14,p(mn)*28,eunit*3,lunit*3,munit*3
     +,pot(mn)*28,grid(mn)*28
      data zero/0.0/


       pi=3.1415 92653 58979d0
      write (*,*)' input prefix of files '
      read (*,99) fname
99    format(a14)
      ln=index(fname,' ')-1
      open(1,file=fname(1:ln)//'.in')
c output file
      open (3,file=fname(1:ln)//'.dm',status='unknown')


      ntypes=0
      nx=0
      ifpair=0
c loop thru reading the input
1     continue
      j=ipickoff(1,p,n,mn)
      if(j.eq.1) go to 2
c repeat input on output
      call echo(3,p,n)


c input energy units and length units
      if(p(1).eq.'UNITS') then
      eunit=p(2)
      lunit=p(3)


      elseif(p(1).eq.'TYPE') then
c input type of particles and hbs2m
       ntypes=ntypes+1
       if(ntypes.gt.2)stop'potgen:ntypes.gt.2'
       hbs2m(ntypes)=rlread(p(3))
c      chrg(ntypes) = rlread(p(4))
 
      elseif(p(1).eq.'GRID')then
c input grid parameters
        do i=1,n
        grid(i)=p(i)
        enddo
        ng=n
        nx=intread(p(2))
        call setgrid(mx,nx,rv,grid(3))


      elseif(p(1).eq.'POT') then
c  input type of potential
        ifpair=ifpair+1
c store away the potential specification
        do i=2,n
         pot(i-1)=p(i)
        enddo
        np=n-1
        if(p(2).eq.'COULOPT'.or.p(2).eq.'SCCLOPT'.
     &     or.p(2).eq.'SCCOULEW') then
         open (4,file=fname(1:ln)//'.vk',status='unknown')
         open (44,file=fname(1:ln)//'.vkk',status='unknown')
        endif


      endif
      go to 1


2     continue
c check if potential and grid have been specified.
      if(ifpair.ne.1) then
       write (*,*) ' ifpair.ne.1'
       stop
      end if
      if(nx.le.0) then
       write (*,*) ' nx.le.0'
       stop
      end if
c zero potential at cutr
      cuto=rlread(pot(2))
      cutr=rlread(pot(3))
       write (*,*) 'real space cut-off:',cutr
       call ipot(pot,np,pot0,cutr,eunit,lunit,ntail,expon,cexpon,
     &           alpha,dummy)
      write (3,33)pot0,ntail,(expon(i),cexpon(i),i=1,ntail)
33    format('POTTAIL ',e14.5,i3,5(f7.2,e13.5))


c setup potential
      do i=1,nx
       call ipot(pot,np,pott(i,1),rv(i),eunit,lunit,ntail,expon,cexpon,
     &           alpha,pott(i,2))
       if(p(2).eq.'MSCCOUL'.or.p(2).eq.'MYUKA') then
        if (rv(i).le.cutr) then
         call ipot(pot,np,vm,2.d0*cutr-rv(i),eunit,lunit,ntail,expon
     &            ,cexpon,alpha,vm2)
         pott(i,1)=pott(i,1)+vm-2.d0*pot0
         pott(i,2)=pott(i,2)+vm2 
        else
         pott(i,1)=0.d0
         pott(i,2)=0.d0
        endif
       endif
      enddo
      write (*,*) 'alpha = ',alpha


c now check grid
      npts=10
      do i=1,nx-1
       diff(i)=0.
       do k=1,npts
         x=rv(i)+(rv(i+1)-rv(i))*k/(npts+1.)
         call ipot(pot,np,potx,x,eunit,lunit,ntail,expon,cexpon,
     &             alpha,dummy)   
         if(p(2).eq.'MSCCOUL'.or.p(2).eq.'MYUKA')then
          if(x.le.cutr) then
           call ipot(pot,np,vm,2.d0*cutr-x,eunit,lunit,ntail,expon
     &              ,cexpon,alpha,dummy)
           potx=potx+vm-2.d0*pot0
          else
           potx=0.d0
          endif
         endif
         call interp(x,ix,a1,a2,a3,a4)
         dp=a1*pott(ix-1,1)+a2*pott(ix,1)+a3*pott(ix+1,1)+
     &      a4*pott(ix+2,1)-potx
         diff(i)=max(diff(i),abs(dp))
       enddo
         diff(nx)=diff(nx-1)
      enddo


      if (p(2).eq.'SCCOUL'.or.p(2).eq.'SCCLOPT'
     &    .or.p(2).eq.'SCCOULEW'.or.p(2).eq.'MSCCOUL'.or.p(2).eq.'MYUKA'
     &   ) then
       rs=rlread(p(5))
       aktf=(12./3.1415926535)**(1./3.)/sqrt(rs)
       write (*,*) 'Thomas-Fermi screening wv',aktf
      endif
      open (11,file=fname(1:ln)//'.pot',status='unknown')
      if (p(2).eq.'SCCOUL'.or.p(2).eq.'SCCOULEW') then
       do i=1,nx
        if (rv(i).le.cutr) then
         write (11,1010)rv(i),pott(i,1)-pot0,pott(i,1)
     &                  ,1./rv(i),exp(-rv(i)*aktf)/rv(i)
        else
         write (11,1010) rv(i),zero,pott(i,1),1./rv(i)
     &                  ,exp(-rv(i)*aktf)/rv(i)
        end if
       end do
      elseif(p(2).eq.'MSCCOUL'.or.p(2).eq.'MYUKA') then
       do i=1,nx
        write (11,1010)rv(i),pott(i,1),1./rv(i),exp(-rv(i)*aktf)/rv(i)
       enddo
      else
       write (*,*) 'v(cutr)=',pot0
       do i=1,nx
        if (rv(i).le.cutr) then
         write (11,1010) rv(i),pott(i,1),pott(i,1)-pot0,pott(i,2)
c        write (11,1010) rv(i),pott(i,1),pott(i,1)-pot0,1./rv(i)
        else
         write (11,1010) rv(i),pott(i,1),zero,pott(i,2)
c        write (11,1010) rv(i),pott(i,1),zero,1./rv(i)
        end if
       end do
      endif
      close (11)
1010  format (5e16.7)
      if(p(2).ne.'MSCCOUL'.and.p(2).ne.'MYUKA') then
       do i=1,nx
        if(rv(i).le.cutr) then
         pott(i,1)=pott(i,1)-pot0
        else
         pott(i,1)=0
         pott(i,2)=0
        endif
       enddo
      endif
      if (p(2).ne.'SCCOULEW')  then
       write(*,*)'do you want a cubic spline at the cut-off? [0=n,1=y]'
       read(*,*) idec
       if (idec.eq.1) then
        call cubspline(nx,rv,pott(1,1),cutr)
c       call cubspline(nx,rv,pott(1,2),cutr)
       endif
      endif
 
cpierleo 16.06.99: FT of the short range SCCOULEW potential
      if (p(2).eq.'SCCOULEW') call epplr(p,alpha,nx,rv,pott(1,1))


      grid(2)='1 '
      write (3,*)'RANK 2 ',nx,1
      call echo(3,grid,ng)
      write (3,*)'LABEL 1 r'
      write (3,*)'BEGIN potential 0'
      write (3,969) (pott(i,1),i=1,nx)
!     grid(2)='1 '
!     write (3,*)'RANK 2 ',nx,1
!     call echo(3,grid,ng)
!     write (3,*)'LABEL 1 r'
!     write (3,*)'BEGIN force 0'
!     write (3,969) (pott(i,2),i=1,nx)
969   format(4e16.7)
      close(3)


      open (9,file=fname(1:ln)//'.dg',status='unknown')
c     write (9,*)'RANK 1 ',nx
c     call echo(9,grid,ng)
c     write (9,*)'LABEL 1 r'
c     write (9,*)'BEGIN grid error'
c     write (9,969) (diff(i),i=1,nx)
      do i=1,nx
      write (9,969) rv(i), diff(i),pott(i,1)
      enddo
      close(9)
cpierleo 12.11.96: sine FFT of the potential
c     cint=0.0
c     do i=2,nx
c      ypot(i)=rv(i)*pott(i,1)
c      cint=cint+rv(i)**2*pott(i,1)
c     enddo
c     cint=cint-.5*rv(nx)**2*pott(nx,1)
c     cint=4.*pi*cint*rv(nx)/nx
c     write (*,*) 'cint=',cint
c     write (12,'(i5,2g20.10)') (i,rv(i),ypot(i),i=1,nx)
c     call sinft(ypot,nx)
c     open (11,file=fname(1:ln)//'.ft',status='unknown')
c     deltan=pi/rv(nx)
c     do i=2,nx
c      rk(i)=deltan*(i-1)
c      ypot(i)=4.*pi*rv(nx)*ypot(i)/rk(i)/nx
c      write (11,1010) rk(i),ypot(i)
c     enddo
c     close(11)
c     hbar=sqrt(hbs2m(1)*8.)
c     do i=1,nx
c      if (rk(20*i).le.40.) then
c       y1=ypot(20*i)
c      else 
c       y1=1.7722e8/rk(20*i)**4
c      endif
c      if (rk(i).le.40.) then
c       y2=ypot(i)
c      else
c       y2=1.7722e8/rk(i)**4
c      endif
c      omega=2./hbar*(1./pi/hbar)**3*
c    &     abs(1.6E5*y1*sin(rk(20*i)*rv(120))-y2*sin(rk(i)*rv(120)))
c      write (13,'(i5,2g20.10)') i,rk(i),omega
c     enddo
      
      stop
      end


      subroutine ipot(p,n,pott,rv,eunit,lunit,ntail,expon,cexpon,
     &                alp,pott2)
      implicit real*8 (a-h,o-z)


      parameter (mnkv=5000,mnsh=900,mdim=3,npts=34)
      integer kmult(0:mnsh),ln,indx
      real*8 tpiell(mdim),rkcomp(mdim,mnkv),rknorm(mnsh),b(40),work(550)
     +  ,vkbare(mnsh),wtk(mnsh),ell(0:mdim),cent,epsil(mnsh)
      real*8 aa,r0,pb,ab(9,4),rr(10)
      real*8 rgrid(npts),vgrid(npts)
      real*8 dr
      parameter (dr=0.1d0)



c interprets potential parameter and sets up table.
      character p(n)*(*),eunit*(*),lunit*(*),cite*80
      dimension expon(10),cexpon(10)
      save icall,nderv,mpoly,cutr,b,c1ry,c2ry,cuto
      save eps,rs,nlamb,rm,a,d,accuracy
      save c6,c8,c9,c10,alpha,beta,gamma,abohr,hartok,cite
      save aa,r0,pb, ab,rr,ln,yuk,hk,x0
      data icall/0/
c helium-helium (aziz-3: MolPhys1987) potential.
c       data eps,rm/10.948,2.963/
c       data a,alpha,beta,d,c6,c8,c10/184431.01,10.43329537,-2.27965105
c    +  ,1.4826,1.36745214,.42123807,.17473318/
c helium-helium (aziz-2) potential
c     data d,eps,a,alpha,c6,c8,c10/1.241314,10.8,544850.4,13.353384
c    +,1.3732412,.4253785,.178100/
c     data rm,beta/2.9673,0./
c helium-helium (aziz-1992) potential.
c      data eps,rm/10.94,2.970/
c      data a,alpha,beta,d,c6,c8,c10/192215.29,10.73520708,-1.89296514
c    +  ,1.4135,1.34920045,.41365922,.17078164/
c exp-six potential for helium from Ross-Young, Phys.Lett.118,463 (1986)
c     data epsry,alphary,rmry/10.8,13.1,2.9673/
c copper many-body potential
      data aa,r0,pb/4.14d-2,2.5565446d0,12.8d0/
      data ab(1,1),ab(2,1),ab(3,1),ab(4,1),ab(5,1),ab(6,1),ab(7,1)
     .,ab(8,1),ab(9,1)/10.652616d0,6.435703d0,0.594270d0,-0.061142d0
     .,-0.054379d0,-0.002970d0,0.015647d0,0.008750d0,0.002267d0/
      data ab(1,2),ab(2,2),ab(3,2),ab(4,2),ab(5,2),ab(6,2),ab(7,2)
     .,ab(8,2),ab(9,2)/-53.263079d0,-33.296982d0,-2.275398d0
     .,-0.357786d0
     .,0.223550d0,0.145520d0,0.018353d0,-0.064040d0,-0.045689d0/
      data ab(1,3),ab(2,3),ab(3,3),ab(4,3),ab(5,3),ab(6,3),ab(7,3)
     .,ab(8,3),ab(9,3)/133.157698d0,66.503273d0,2.433580d0,1.045397d0
     .,0.097167d0,-0.409287d0,-0.017448d0,-0.130167d0,0.306888d0/
      data ab(1,4),ab(2,4),ab(3,4),ab(4,4),ab(5,4),ab(6,4),ab(7,4)
     .,ab(8,4),ab(9,4)/-222.181416d0,-47.459032d0,-0.839491d0
     .,-0.621220d0
     .,-0.675271d0,0.438298d0,-0.067316d0,1.402978d0,-0.687109d0/
      data rr(1),rr(2),rr(3),rr(4),rr(5),rr(6),rr(7),rr(8),rr(9)
     ./1.4500d0,1.550d0,2.d0,2.5512d0,3.06d0,3.31d0,3.608d0,4.1662d0
     .,4.27d0/
! Kolos-Wolniesky potential for H2
      data (rgrid(j),j=1,6)/0.40, 0.50, 0.60, 0.70, 0.80, 0.90/
      data (rgrid(j),j=7,12)/1.00, 1.10, 1.20, 1.30, 1.40, 1.50/
      data (rgrid(j),j=13,18)/1.60, 1.70, 1.80, 1.90, 2.00, 2.10/
      data (rgrid(j),j=19,24)/2.20, 2.30, 2.40, 2.50, 2.60, 2.70/
      data (rgrid(j),j=25,30)/2.80, 2.90, 3.00, 3.10, 3.20, 3.30/
      data (rgrid(j),j=31,34)/3.40, 3.50, 3.60, 3.70/
      data (vgrid(j),j=1,3)/-0.1202028, -0.5266270, -0.7696253/
      data (vgrid(j),j=4,6)/-0.9220185, -1.0200487, -1.0836362/
      data (vgrid(j),j=7,9)/-1.1245331, -1.1500512, -1.1649294/
      data (vgrid(j),j=10,12)/-1.1723414, -1.1744699, -1.1728492/
      data (vgrid(j),j=13,15)/-1.1685773, -1.1624521, -1.1550616/
      data (vgrid(j),j=16,18)/-1.1468425, -1.1381236, -1.1291528/
      data (vgrid(j),j=19,21)/-1.1201190, -1.1111659, -1.1024035/
      data (vgrid(j),j=22,24)/-1.0939149, -1.0857627, -1.0779927/
      data (vgrid(j),j=25,27)/-1.0706404, -1.0637259, -1.0572607/
      data (vgrid(j),j=28,30)/-1.0512547, -1.0457057, -1.0406020/
      data (vgrid(j),j=31,33)/-1.0359419, -1.0318402, -1.0278471/
      data (vgrid(j),j=34,34)/-1.0243742/



      icall=icall+1
      pott2=0.0
      if(icall.eq.1) then
       ntail=0
!      ln=index(p(1),' ')-1
!      write(*,'(a)') p(1)(1:ln)
      endif


      if(p(1).eq.'HEDF2') then
       alp=0.0
       if(icall.eq.1) then
c       xe=convert('energy','meV',eunit)
c       xl=convert('length','au',lunit)
c He-He potential Aziz, PRL 74, 1586 (1995) ab initio HFD-B3-FCI1
        eps=10.956d0
        rm= 2.9683d0
        a=1.86924404d5
        alpha=10.5717543d0
        beta=-2.07758779d0
        d=1.438d0
        c6=1.35186623d0
        c8=.41495143d0
        c10=.17151143d0
        cite='R. A. Aziz et al. , PRL 74,1586 (1995)'
        ntail=3
        expon(1)=-6.
        cexpon(1)=-c6*eps*rm**6
        expon(2)=-8.
        cexpon(2)=-c8*eps*rm**8
        expon(3)=-10.
        cexpon(3)=-c10*eps*rm**10
        write (3,33)cite
33      format('CITATION:',a80)
        
       endif
 
        x=rv/rm
        a2i=1./(x*x)
        a6i=a2i*a2i*a2i
        if(x.lt.d) then
          f=exp(-(d/x-1.)**2)
          f2=2.*(d/x-1.)*d*exp(-(d/x-1.)**2)/x**2
        else
          f=1.d0
          f2=0.d0
        endif
         pott=eps*(a*exp(-alpha*x+beta*x*x)-f*a6i*
     +  (c6+a2i*(c8+a2i*c10)))
         pott2=-eps*(a*(-alpha+2*beta*x)*exp(-alpha*x+beta*x*x)-
     &              f2*a6i*(c6+a2i*(c8+a2i*c10))+
     &              6.*f*a6i*(c6+a2i*(c8+a2i*c10))/x+
     &              2.*f*a6i*(a2i*(c8+2.*a2i*c10))/x)/rm


      elseif(p(1).eq.'HE86RY') then
       if (icall.eq.1) then
        c1ry=epsry*6./(alphary-6.)
        c2ry=-epsry*alphary/(alphary-6.)
       endif
       x=rmry/rv
       x1=1./x
       x6=x**6
       pott=c1ry*exp(alphary*(1.-x1))+c2ry*x6


      elseif(p(1).eq.'HeNe') then
       echg=1.6021917e-19
       dkbolt=1.380622e-23
       evtokb=echg/dkbolt
      rm=3.029
      sig=2.699
      eps=0.001827*evtokb
      alfa=4.031
      beta=0.0987
      a=797.*evtokb
      d=1.28
      c6=1.810*evtokb
      c8=5.503*evtokb
      c10=20.5*evtokb
      c12=96.*evtokb
      c14=554.*evtokb
      drm=d*rm
       if(icall.eq.1) then
        ntail=5
        expon(1)=-6.
        cexpon(1)=-c6
        expon(2)=-8.
        cexpon(2)=-c8
        expon(3)=-10.
        cexpon(3)=-c10
        expon(4)=-12.
        cexpon(4)=-c12
        expon(5)=-14.
        cexpon(5)=-c14
       endif
 
      r=rv
      vhf=a*exp(-r*(alfa+r*beta))
      ri2=1./(r*r)
      ri6=ri2*ri2*ri2
      vd=-ri6*(c6+ri2*(c8+ri2*(c10+ri2*(c12+ri2*c14))))
      if(r.lt.drm)then
        f=drm/r-1.
        vd=vd*exp(-f*f)
      endif
      pott=vhf+vd


      elseif(p(1).eq.'WALL') then
        rc=rlread(p(2))
        if(n.ge.3) then
          cent=rlread(p(3))
        else
          cent=0.
        endif
        pottmin=1000.
        ntail=0
        x=(rc-rv)/rm
         if(x.le.0)then
          pott=pottmin
         else
        a2i=1./(x*x)
        a6i=a2i*a2i*a2i
        if(x.lt.d) then
          f=exp(-(d/x-1.)**2)
        else
          f=1
        endif
         pott=eps*(a*exp(-alpha*x+beta*x*x)-f*a6i*
     +  (c6+a2i*(c8+a2i*c10)))
         if(pott.lt.0)pott=0.
         if(pott.gt.pottmin)pott=pottmin
         endif
         pott=pott+cent/rv**2
         
 
      elseif(p(1).eq.'COUL') then
       nexp=1
       cutr=rlread(p(2))
       eps=rlread(p(3))
       accuracy=rlread(p(4))


       if(accuracy.gt.0) then
        alpha=-.5*log(cutr*accuracy/abs(eps))/cutr**2
        xpo=.5*nexp
        call gammi(gr,xpo,alpha*rv**2,gr0)
        pott=eps*gr/(gr0*rv**nexp)
        alp = alpha
       else
        pott=eps/rv**nexp
       endif


      elseif(p(1).eq.'MYUKA') then
       cutr=rlread(p(2))
       eps=rlread(p(3))
       tfk=rlread(p(4))
       pott=eps*dexp(-tfk*rv)/rv


      elseif(p(1).eq.'SCCOUL'.or.p(1).eq.'MSCCOUL') then
       if(icall.eq.1)then
        eps=rlread(p(3))
        rs=rlread(p(4))
        alp=0.0
       endif
       call epp(eps,rs,rv,pott)


      elseif(p(1).eq.'SCCOULEW') then
       if(icall.eq.1)then
        cutr=rlread(p(2))
        eps=rlread(p(3))
        rs=rlread(p(4))
        accuracy=rlread(p(5))
        alpha=-.5*log(cutr*accuracy/abs(eps))/cutr**2
        alp=alpha
       endif   
       call epp(eps,rs,rv,pott)
       call gammi(gr,.5d0,alpha*rv**2,gr0)
       pott=pott*eps*gr/gr0


c new option for screened-optimized breakup coulomb
      elseif(p(1).eq.'SCCLOPT') then
       if(icall.eq.1)then
         if(n.lt.5) then
          write (*,*)' two few parameters for SCCLOPT',n
          stop
         endif
       cutr=rlread(p(2))
       eps=rlread(p(3))
       rs=rlread(p(4))
       ndim=3  ! three dimensions only
       nlamb=intread(p(5))
       pi=3.1415 92653 58979d0
c default with 5 paramters is cubic box and ell=2*cutr. Otherwise read in ell.
       vol=1.
       do l=1,ndim
        if(l+5.le.n) then
         ell(l)=rlread(p(l+5))
        else
         ell(l)=ell(l-1)
        endif
        vol=vol*ell(l)
        if(1.9999*cutr.gt.ell(l))stop 'ipot: 1.9999*cutr.gt.ell(l)'
        tpiell(l)=2*pi/ell(l)
       enddo
c      do l=1,ndim
c       ell(l)=2.d0*cutr
c       tpiell(l)=2*pi/ell(l)
c      enddo
c       !NOTE WE ASSUME ELL=2*cutr
c      vol= ell(1)**ndim
c this is the order of the polynomial
!        mpoly=9
c number of derviatives zeroed at cutr
!        nderv=2
c cutoff of explicit k summation
         cutk = 6*tpiell(1)
         cutk=18.849556d0/cutr
c cutoff for approximate k summation
777      argek=20*cutk
         write (*,*) 'small k=',cutk,' large k=',argek,' mpoly=',mpoly,
     +    ' nderv=',nderv
c generate the set of kvectors
         call shells(ndim,tpiell,cutk,nshlls,rkcomp,rknorm(2),kmult
     +      ,nvects,mnkv,mnsh,mdim)
         if (nshlls.lt.nlamb) then
          cutk=1.1*cutk
          goto 777
         endif
         call fillk(rknorm,wtk,kmult,nshlls,nshex1,argek,ndim,vol,mnsh)
         sangle=sang(ndim)
        !get epsilons for screening
cpierleo   WARNING:
c the routine epsk1 works in units of a therefore the reciprocal 
c distances must be converted before going into that.
c
         if (lunit.ne.'A') then
          if (lunit.eq.'A0') then
           do k=1,nshex1
            rknorm(k) = rknorm(k)/rs
           enddo
          else
           write (*,*) 'unit of length not supported: lunit=',lunit
           stop
          endif
         endif
         call epsk1(rs,epsil(2),rknorm(2),nshex1,gm) ! evaluate epsilon
         if (lunit.ne.'A') then
          if (lunit.eq.'A0') then
           do k=1,nshex1
            rknorm(k) = rknorm(k)*rs
           enddo
          else
           write (*,*) 'unit of length not supported: lunit=',lunit
          endif
         endif
!        write (99,'(2g18.8)') (rknorm(ik),epsil(ik),ik=1,nshex1)
        !this is fourier transform of potential
         do k=2,nshex1
          vkbare(k)=(eps*sangle)/(vol*rknorm(k)**(ndim-1)*epsil(k))
         enddo
         write (*,*)'gamma0=',gm
c  vkbare(1)=[(k_f/k_tf)^2/rs-gamma0]/k_f^2 is the k=0 limit of the potential
         vkbare(1)=((3*pi**2/16)**(2./3.)/rs-gm)/(1.919158293/rs)**2
         vkbare(1)=(eps*sangle)*vkbare(1)
         write (*,*) 'limiting value of the potential for small k:'
     &                ,vkbare(1)
         vkbare(1)=vkbare(1)/vol
c give the interaction energy with images:
        vmad=vzero(rs)
        write (*,*) ' vmad',vmad

        !determine best polynomial
!        call fitpn(mpoly,vkbare,rknorm,wtk,nshex1,nlamb,eps
!    +      ,b,work,0,vol,cutr,ndim,1,nderv,vmad)
!FITPNNEW CODING
          nknots=30
          m=2
          maxn=2*m+1
          vt0=eps
          vt1=0.d0
          coul=.true.
          t0=.true.   ! fixes const term for r=0 with value vt0
          t1=.false.
         write (*,*) 'small k=',cutk,' large k=',rknorm(nshex)
     & ,' knots=',nknots,' nderv=',nderv
          call fitpnnew(vkbare,vlr,b,rknorm,wtk,nshex,nlamb,m,maxm
     &   ,cutr,nknots,t0,vt0,t1,vt1,coul,vmad,vol,sm,rknot,delta,ndim)
         write (*,*)' finished fitpn'


c Vimage is 1/2 limit (r=0) [ V_image(r) -V_bare(r) ]
         write (3,*)'VIMAGE ',.5*vmad
c write out k-space potential. Unit 4 should have been opened already.
          write (4,*)'RANK 1 ',nlamb+1
          write (4,*)'BEGIN k-space potential'
!         write (4,*) (vkbare(k),k=1,nlamb+1)
          write (4,*) (vlr(k),k=0,nlamb)
          close(4)
!         write (44,'(15H#  k     Yk     )')
!         write (44,'(3g20.10)') 
!    &          (rknorm(k),vkbare(k),wtk(k),k=2,nlamb+1)
!         close(44)
     
       endif
c now evaluate potential
      call computespl2(rv,pott,dpot,ddpot,d3pot,m
     .           ,maxm,maxn,nknots,sm,rknot,b,delta)
      pott=pott/rv
!      x=rv/cutr
!      pott=0.
!      do i=mpoly-nderv,1,-1
!          pott=b(i)+x*pott
!      enddo
!      pott=pott*(x-1.)**nderv/x
! end of scclopt
      elseif(p(1).eq.'COULOPT') then
       if(icall.eq.1) then
         if(n.lt.5) then
          write (*,*)' two few parameters for COULOPT',n
          stop
         endif
         cutr=rlread(p(2))
         eps=rlread(p(3))
         ndim=intread(p(4))
         if(ndim.lt.2.or.ndim.gt.3)stop
         nlamb=intread(p(5))
         if(n.eq.5) ell(0)=2*cutr
         vol=1.
         pi=3.1415 92653 58979
c default with 5 paramters is cubic box and ell=2*cutr. Otherwise read in ell.
         do l=1,ndim
           if(l+5.le.n) then
            ell(l)=rlread(p(l+5))
           else
            ell(l)=ell(l-1)
           endif
           vol=vol*ell(l)
           if(1.9999*cutr.gt.ell(l))stop
           tpiell(l)=2*pi/ell(l)
         enddo
 
c cutoff of explicit k summation
         cutk = 6*tpiell(1)
         cutk=13.d0/cutr 
c cutoff for approximate k summation
         argek=50*cutk
c this is the order of the polynomial
         mpoly=8
c number of derviatives zeroed at cutr
         nderv=2
         write (*,*) 'small k=',cutk,' large k=',argek,' mpoly=',mpoly,
     +    ' nderv=',nderv
c generate the set of kvectors
         call shells(ndim,tpiell,cutk,nshlls,rkcomp,rknorm(2),kmult
     +      ,nvects,mnkv,mnsh,mdim)
         call fillk(rknorm,wtk,kmult,nshlls,nshex1,argek,ndim,vol,mnsh)
c put in uniform background by setting k=0 term to zero
         sangle=sang(ndim)
         vkbare(1)=0.
         do k=2,nshex1
            vkbare(k)=(eps*sangle)/(vol*rknorm(k)**(ndim-1))
         enddo 
 
c Compute madeulung constant.
          xpon=1
         call mdlng(ndim,xpon,ell(1),eps,vmad)
c Find optimal division into k and r space.
       write (*,*)' vmad ',vmad


       if(ndim.eq.3)vmad=-2.837297479*eps/ell(1)
       if(ndim.eq.2)vmad=-3.90026492*eps/ell(1)
       write (*,*) ' ideal square madelung constant',vmad


          do k=1,nshex1
           write(94,*)k,wtk(k),rknorm(k),vkbare(k)
          enddo
            write(94,*)nshex1,nlamb,eps,0,vol,cutr,ndim,1
     .                ,nderv,vmad
         call fitpn(mpoly,vkbare,rknorm,wtk,nshex1,nlamb,eps
     +      ,b,work,0,vol,cutr,ndim,1,nderv,vmad)
c Vimage is 1/2 limit (r=0) [ V_image(r) -V_bare(r) ]
         write (3,*)'VIMAGE ',.5*vmad
c write out k-space potential. Unit 4 should have been opened already.
          write (4,*)'RANK 1 ',nlamb+1
          write (4,*)'BEGIN k-space potential'
          write (4,*) (vkbare(k),k=1,nlamb+1)
          close(4)
          write(94,*) 'dopo fitpn'
          do k=1,nlamb+1
           write(94,*) rknorm(k),vkbare(k)
          enddo
          close(94)
       endif


c now evaluate potential
       x=rv/cutr
       pott=0.
       do i=mpoly-nderv,1,-1
           pott=b(i)+x*pott
       enddo
       pott=pott*(x-1.)**nderv/x


      elseif(p(1).eq.'LJ') then
       if(icall.eq.1)then
        eps=rlread(p(3))
        sigma=rlread(p(4))
        ntail=2
        expon(1)=-6.
        cexpon(1)=-4*eps*sigma**6
        expon(2)=-12.
        cexpon(2)=4*eps*sigma**12
       endif
        x=rv/sigma
        a2i=1./(x*x)
        a6i=a2i*a2i*a2i
        pott=4*eps*a6i*(a6i-1.)
      elseif(p(1).eq.'POWER') then
       if(icall.eq.1)then
        eps=rlread(p(3))
        sigma=rlread(p(4))
        expon(1)=-rlread(p(5))
        ntail=1
        cexpon(1)=4*eps*sigma**(-expon(1))
       endif
        x=rv/sigma
        pott=4*eps*x**expon(1)


      elseif(p(1).eq.'H2H2') then
       if(icall.eq.1) then
c potential parameters in atomic units
        abohr=0.529177249e00
c conversion factor for energy from hartree to K
        hartok=3.15777321e05
        alpha= 1.713e00
        beta= 1.5671e00
        gamma=0.00993e00
        c6=12.14e00
        c8=215.2e00
        c10=4813.9e00
        c9=143.1e00
        cite=
     +  'I. F. Silvera and V. V. Goldman, J. Chem. Phys. 69, 4209(1978)'
        ntail=4
        expon(1)=-6.
        cexpon(1)=-c6*eps*abohr**6
        expon(2)=-8.
        cexpon(2)=-c8*eps*abohr**8
        expon(3)=-9.
        cexpon(3)=+c9*eps*abohr**9
        expon(4)=-10.
        cexpon(4)=-c10*eps*abohr**10
        write (3,33)cite
       endif


c distance in abohr
       rmm=3.41/abohr
       x=rv / abohr


       ri=1.e00 / x
       if(x .lt. 1.28*rmm) then
        f=exp(-(rmm*1.28 * ri - 1.e00)**2)
       else
        f=1.e00
       endif
       r2i=ri * ri
       r6i=r2i * r2i * r2i


       pott=hartok*( exp(alpha - x * (beta + gamma * x))
     &           - f*r6i*(c6 + r2i*(c8 + r2i*c10 - c9*ri)) )



      elseif (p(1).eq.'HARM') then
       if (icall.eq.1) then
        alpha=rlread(p(2))
       endif
       pott=alpha*rv**2


      elseif(p(1).eq.'COPPER') then
       if (icall.eq.1) then
        cutr=rlread(p(2))
        rr(10)=cutr
       endif
       if(rv.lt.cutr) then
        if (rv.ge.rr(1)) then
         do i=1,9
          if(rv.ge.rr(i).and.rv.le.rr(i+1))then
           pott=0.d0
           do j=1,4
            pott=pott+ab(i,j)*(rv-rr(i))**(j-1)
           enddo
          endif
         enddo
        else
         pott=aa*exp(-pb*(rv/r0-1.d0))
        endif
       else
        pott=0.d0
       endif


      elseif(p(1).eq.'KW') then  ! Kolos-Wolniesky H2 potential
       if(icall.eq.1) then
        cuto=rlread(p(2))
        cutr=rlread(p(3))
        write(93,'(2g20.10)') (rgrid(j),vgrid(j),j=1,34)
       endif
       if (rv.ge.cuto.and.rv.le.cutr) then
        indx = int((rv-rgrid(1))/dr) + 1
        pott = (vgrid(indx+1)-vgrid(indx))*(rv-rgrid(indx))/dr
     &            +  vgrid(indx)
       else
        pott=0.d0
       endif
 
      elseif(p(1).eq.'HBOND') then  
       if(icall.eq.1) then
        yuk=rlread(p(2))
        cutr=rlread(p(3))
        hk=rlread(p(4))
        x0=rlread(p(5))
       endif
       if(rv.le.cutr)then
        pott=dexp(-yuk*rv)/rv+dexp(-yuk*(2.d0*cutr-rv))/(2.d0*cutr-rv)
     &      +0.5d0*hk*(rv-x0)**2
       else
        pott=0.d0
       endif


      else
         write (*,*)' this potential not defined ',p(1),n
         stop


      endif


      return
      end


      subroutine epp(q,rs,r,phir)
      implicit real*8 (a-h,o-z)
c finds the interaction between 2 protons in an electron gas
c rs is the electon gas density, the potential phir will be computed
c on the grid r. r is in atomic units. k in units of rs.
      parameter (nk=10000,dk=.01)
      dimension phik(nk),eps(nk)
      save pi,icall,phik,rsi
      data pi,icall/3.1415926535,0/


      if(icall.eq.0) then
c calculate fourier transform of phi(r)*r first time through
       call epsk(rs,eps,dk,nk)
       akn=dk/(.75*pi)
       rsi=1./rs


       do k=1,nk
        phik(k)=akn/(eps(k)*k*dk)
        if(mod(k,2).eq.1) phik(k)=2*phik(k)
       enddo
       phik(nk)=.5*phik(nk)
      endif


      icall=icall+1


      z=nk*dk*r*rsi
      zi=1./z
c  corrections for the FT at large k.
      phir=(2./pi)*zi*(cos(z)*(1.-2.*zi**2)+sin(z)*zi)


      do k=1,nk
      rk=k*dk*rsi
      phir=phir+phik(k)*sin(rk*r)
      enddo
 
      if(r*dk*nk.lt.2.) phir=1.
      phir=q*phir/r
      return
      end


      subroutine epsk(rs,eps,dk,nk)
      implicit real*8 (a-h,o-z)
      save pi,fwv,gm
      dimension eps(nk)
      data pi/3.1415926535/
      data fwv /1.919158293/
      alp=.082929555*rs
c alp=[(k_TF/k_F)**2]/8
      do 1 k=1,nk
      rk=k*dk
      eta=rk/(2*fwv)
      ff=1+(1-eta**2)*log((1+eta)/abs(1-eta))/(2.*eta)
1     eps(k)=1+1/(-g(rk,rs,ff,gm)+eta**2/(alp*ff))
      return
      end


      subroutine epsk1(rs,eps,rk,nk,gm)
      implicit real*8 (a-h,o-z)
      save pi,fwv
      dimension eps(nk),rk(nk)
      data pi/3.1415926535d0/
      data fwv /1.919158293d0/  ! k_F/a
      alp=.082929555d0*rs       ! 1/8*(k_TF/k_F)^2*rs
c
c  -v(k)*Chi_L(k)=alp*ff/eta^2
c
      do 1 k=1,nk
      eta=rk(k)/(2.d0*fwv)
      ff=1+(1-eta**2)*log((1+eta)/abs(1-eta))/(2.d0*eta)
1     eps(k)=1.d0+1.d0/(-g(rk(k),rs,ff,gm)+eta**2/(alp*ff))
      return
      end


      function vzero(rs)
      implicit real*8 (a-h,o-z)
      parameter (nk=20000,dk=.001)
      save pi,fwv
      data pi/3.1415926535d0/
      data fwv /1.919158293d0/  ! k_F/a
      alp=.082929555d0*rs       ! 1/8*(k_TF/k_F)^2*rs
c
c  -v(k)*Chi_L(k)=alp*ff/eta^2
c
      vzero=.5
      do  k=1,nk
       rk=k*dk
       eta=rk/(2.d0*fwv)
       ff=1+(1-eta**2)*log((1+eta)/abs(1-eta))/(2.d0*eta)
       gg=-g(rk,rs,ff,gm)+eta**2/(alp*ff)
       gg=1./(gg+1.)
       vzero=vzero+gg
       if (mod(k,100).eq.0) write (66,*) rk,gg
       enddo
       ! assuming that gg drops off as rk^{-4}
      vzero=(2./pi)*(dk*vzero+gg*(-.5*dk+rk/3.))
      return
      end


 
      function ecvwm(rs)
      implicit real*8 (a-h,o-z)
c correlation energy of the electron gas
c from relation [4.4] of ref: Vosko,Wilk and Nusair, Can. J. P. 58, 1200 (1980)
      save x0,b,c,a
      data x0,b,c/-.10498,3.72744,12.9352/
      data a/.0621814/
      xf(x)=c+x*(b+x)
      q=sqrt(4*c-b**2)
      x=sqrt(rs)
      xc=xf(x)
      xc0=xf(x0)
      t=2*atan(q/(2*x+b))/q
      ecvwm=a*(log(x**2/xc)+b*t-b*x0*(log((x-x0)**2/xc)+(b+2*x0)*t
     +)/xc0)
      return
      end


      function ehf(rs)
      implicit real*8 (a-h,o-z)
      ehf=2.2099/rs**2-.91633/rs
      return
      end


  
      function g(rk,rs,ff,gm)
      implicit real*8 (a-h,o-z)
      save a,b,c,b1,rls,fwv,fwvi,gmm
c local field correction from Ichimaru, Iyetomi and Tanaka, Phys. Rep. 149, 
c                                                           91 (1987).
c  see page 121-122 for the explanation of this routine
      data rls/0./
      data fwv /1.919158293/   ! this is the fermi weave vector in angstrom^-1
      if(rls.ne.rs) then
c first time through find a b c b1
      ec=ecvwm(rs)
      ecp=ecvwm(rs+1.e-5)
      ecm=ecvwm(rs-1.e-5)
      de=.5e5 * (ecp-ecm)
      dde= 1.e10 * ( (ecp-ec) + (ecm-ec) )
      gmm=.25-.0682068*rs**2*(rs*dde-2*de)
      g0=.125/bi1(1.62903*sqrt(rs))**2
      a=.029
      b=.5625*gmm-3.*(1.-g0)/64. -16.*a/15.
      c=-.75*gmm+.5625*(1.-g0)-3.2*a
      b1=b+8*a/3.
      fwvi=1./fwv
      write (*,*) 'local field correction parameters:'
      write (*,'(13hgamma0, g(0):,2g18.8)') gmm,g0
      write (*,'(8hA, B, C:,3g18.8)') a, b, c     
      rls=rs
      endif


      q=rk*fwvi
      q2=q**2
      g=c+q2*(b+a*q2)+(-c+q2*(b1+q2*a))*(ff-1.)
      gm=gmm
      return
      end
 
      subroutine epplr(p,alpha,nx,rv,v)
      implicit none
      integer nk,nx,i,j,k
      real*8 pi
      parameter(nk=250,pi=3.1415 92653 58979d0)
      real*8 v(nx),rv(nx),eps,rs,dk,fpi,vsk(nk),epsil(nk)
      real*8 vkl(nk),rlread,rkmax,rk(nk),vmad,vzero,alpha
      real*8 vkbare(nk),gm,dr
      character p(10)*28
      save gm 
      
      eps=rlread(p(4))
      rs=rlread(p(5))
      fpi=4.d0*pi
      dr=rv(2)-rv(1)
      rkmax=0.5d0/dr
      dk=rkmax/nk
c     do i=1,nx
c      v(i)=dexp(-rv(i)**2)
c     enddo
      rk(1)=0.d0
      do k=2,nk
       rk(k)=(k-1)*dk
       vsk(k)=rk(k)*(1.d0-dcos(rk(k)*rv(1)))
c      vsk(k)=dsin(rk(k)*rv(1))/rk(k)-rv(1)*dcos(rk(k)*rv(1))
       vsk(k)=vsk(k)+.5d0*v(1)*rv(1)*dsin(rv(1)*rk(k))*dr
       do i=2,nx-1
        vsk(k)=vsk(k)+v(i)*rv(i)*dsin(rv(i)*rk(k))*dr
       enddo
       vsk(k)=vsk(k)+0.5d0*v(nx)*rv(nx)*dsin(rv(nx)*rk(k))*dr
       vsk(k)=fpi*eps*vsk(k)/rk(k)
      enddo 
      do k=1,nk
       rk(k)=rk(k)*rs
      enddo
      call epsk1(rs,epsil(2),rk(2),nk-1,gm)
      do k=1,nk
       rk(k)=rk(k)/rs
      enddo
      do k=2,nk
       vkbare(k)=fpi*eps/rk(k)**2/epsil(k)
      enddo
c vsk(1) is the k=0 fourier transform of the short range pot.
      vsk(1)=rv(1)**2*(1.d0+v(1)*dr)/2.d0
      do i=2,nx-1
       vsk(1)=vsk(1)+rv(i)**2*v(i)*dr
      enddo
      vsk(1)=vsk(1)+v(nx)*rv(nx)**2/2.d0*dr
      vsk(1)=vsk(1)*fpi*eps
      write(*,*) 'gamma0=',gm
c  vkl(1)=[(k_f/k_tf)^2/rs-gamma0]/k_f^2-vsk(0) 
c  is the k=0 limit of the long range potential
      vkbare(1)=((3*pi**2/16)**(2./3.)/rs-gm)/(1.919158293/rs)**2
      vkbare(1)=fpi*eps*vkbare(1)
      write (*,*) 'limiting value of the potential at k=0:'
     &                ,vkbare(1)-vsk(1)
      do k=1,nk
       vkl(k)=vkbare(k)-vsk(k)
      enddo
      write(95,'(4g20.10)') (rk(k),vkbare(k),vsk(k),vkl(k),k=1,nk)
c give the interaction energy with images:
      vmad=-2.d0*dsqrt(alpha/pi)
      write (*,*) ' vmad',vmad
c Vimage is 1/2 limit (r=0) [ V_image(r) -V_bare(r) ]
      write (3,*)'VIMAGE ',.5*vmad
c write out k-space potential. Unit 4 should have been opened already.
      write (4,*)'RANK 1 ',nk
      write (4,*)'BEGIN k-space potential'
      write (4,*) (vkl(k),k=1,nk)
      close(4)
c     write (44,*) '#  k     Vlr(k)  '
      write (44,'(15H#  k     Yk     )')
      write (44,'(2g20.10)') (rk(k),vkl(k),k=1,nk)
      close(44)
      return
      end
