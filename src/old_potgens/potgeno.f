      program potgen
      implicit real*8 (a-h,o-z)
c generates a potential table Updated to include optimized coulomb 4/30/93
      parameter (mn=10,mx=1000)
      dimension pott(mx),rv(mx),diff(mx),expon(20),cexpon(20)
      character fname*14,p(mn)*28,eunit*3,lunit*3
     +,pot(mn)*28,grid(mn)*28

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
        if(p(2).eq.'COULOPT')open (4,file=fname(1:ln)//'.vk'
     +     ,status='unknown')

      endif
      go to 1

2     continue
c check if potential and grid have been specified.
      if(ifpair.ne.1)stop
      if(nx.le.0)stop
c setup potential
      do i=1,nx
      call ipot(pot,np,pott(i),rv(i),eunit,lunit,ntail,expon,cexpon)
      enddo

c now check grid
      npts=10
      do i=1,nx-1
       diff(i)=0.
       do k=1,npts
         x=rv(i)+(rv(i+1)-rv(i))*k/(npts+1.d0)
         call ipot(pot,np,potx,x,eunit,lunit,ntail,expon,cexpon)   
         call interp(x,ix,a1,a2,a3,a4)
         dp=a1*pott(ix-1)+a2*pott(ix)+a3*pott(ix+1)+a4*pott(ix+2)-potx
         diff(i)=max(diff(i),abs(dp))
       enddo
         diff(nx)=diff(nx-1)
      enddo

c zero potential at cutr
      cutr=rlread(pot(2))
       call ipot(pot,np,pot0,cutr,eunit,lunit,ntail,expon,cexpon)
      write (3,33)pot0,ntail,(expon(i),cexpon(i),i=1,ntail)
33    format(' POTTAIL ',e14.5,i3,5(f7.2,e13.5))
       do i=1,nx
       if(rv(i).le.cutr) then
         pott(i)=pott(i)-pot0
       else
         pott(i)=0
       endif
      enddo
 
      grid(2)='1 '
      write (3,*)'RANK 2 ',nx,1
      call echo(3,grid,ng)
      write (3,*)'LABEL 1 r'
      write (3,*)'BEGIN potential 0'
      write (3,969) (pott(i),i=1,nx)
969   format(5e17.9)
      close(3)

      open (9,file=fname(1:ln)//'.dg',status='unknown')
c     write (9,*)'RANK 1 ',nx
c     call echo(9,grid,ng)
c     write (9,*)'LABEL 1 r'
c     write (9,*)'BEGIN grid error'
c     write (9,969) (diff(i),i=1,nx)
      do i=1,nx
      write (9,969) rv(i), diff(i),pott(i)
      enddo
      close(9)
      end

      subroutine ipot(p,n,pott,rv,eunit,lunit,ntail,expon,cexpon)
      implicit real*8 (a-h,o-z)

      parameter (mnkv=1000,mnsh=900,mdim=3)
      integer kmult(0:mnsh)
      real*8 tpiell(mdim),rkcomp(mdim,mnkv),rknorm(mnsh),b(40),work(550)
     +  ,vkbare(mnsh),wtk(mnsh),ell(0:mdim),cent

c interprets potential parameter and sets up table.
      character p(n)*(*),eunit*(*),lunit*(*),cite*80
      dimension expon(10),cexpon(10)
      save icall,nderv,mpoly,cutr,b
      save c6,c8,c9,c10,alpha,beta,gamma,abohr,hartok
      data icall/0/
      icall=icall+1
      if(icall.eq.1)ntail=0

      if(p(1).eq.'HEDF2') then
c helium-helium (aziz-3) potential.
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
c He-He potential Aziz, PRL 74, 1586 (1995) ab initio HFD-B3-FCI1
       data eps,rm/10.956d0,2.9683d0/
       data a,alpha,beta,d,c6,c8,c10/1.86924404d5,10.5717543d0
     +,-2.07758779d0
     +,1.438d0,1.35186623d0,.41495143d0,.17151143d0/
       cite=
     +'R. A. Aziz, A. R. Janzen, M. R. Moldover, Phys. Rev. Letts. 74,
     +1586 (1995)'
       if(icall.eq.1) then
c      xe=convert('energy','meV',eunit)
c      xl=convert('length','au',lunit)
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
        else
          f=1
        endif
         pott=eps*(a*exp(-alpha*x+beta*x*x)-f*a6i*
     +  (c6+a2i*(c8+a2i*c10)))

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
      accuracy=rlread(p(4))
      eps=rlread(p(3))

      if(accuracy.gt.0) then
        alpha=-.5*log(cutr*accuracy/abs(eps))/cutr**2
        xpo=.5d0*nexp
        call gammi(gr,xpo,alpha*rv**2,gr0)
        pott=eps*gr/(gr0*rv**nexp)
       else
        pott=eps/rv**nexp
       endif

      elseif(p(1).eq.'SCCOUL') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       rs=rlread(p(4))
      endif
       call epp(eps,rs,rv,pott)

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
c cutoff for approximate k summation
         argek=10*cutk
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

         call fitpn(mpoly,vkbare,rknorm,wtk,nshex1,nlamb,eps
     +      ,b,work,0,vol,cutr,ndim,1,nderv,vmad)
c Vimage is 1/2 limit (r=0) [ V_image(r) -V_bare(r) ]
         write (3,*)'VIMAGE ',.5*vmad
c write out k-space potential. Unit 4 should have been opened already.
          write (4,*)'RANK 1 ',nlamb+1
          write (4,*)'BEGIN k-space potential'
          write (4,*) (vkbare(k),k=1,nlamb+1)
          close(4)
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
      cite='I. F. Silvera and V. V. Goldman, J. Chem. Phys. 69, 42
     +09(1978)'
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
      parameter (nk=3200,dk=.02)
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
      save pi,fwv
      dimension eps(nk)
      data pi/3.1415926535/
      data fwv /1.919158293/
      alp=.082929555*rs
      do 1 k=1,nk
      rk=k*dk
      eta=rk/(2*fwv)
      ff=1+(1-eta**2)*log((1+eta)/abs(1-eta))/(2.*eta)
1     eps(k)=1+1/(-g(rk,rs,ff)+eta**2/(alp*ff))
      return
      end

 
      function ecvwm(rs)
      implicit real*8 (a-h,o-z)
c correlation energy of the electron gas
      save x0,b,c,a
      xf(x)=c+x*(b+x)
      data x0,b,c/-.10498,3.72744,12.9352/
      data a/.0621814/
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

  
      function g(rk,rs,ff)
      implicit real*8 (a-h,o-z)
      save a,b,c,b1,rls,fwv,fwvi
c local field correction
      data rls/0./
      data fwv /1.919158293/

      if(rls.ne.rs) then
c first time through find a b c b1
      ec=ecvwm(rs)
      ecp=ecvwm(rs+1.e-5)
      ecm=ecvwm(rs-1.e-5)
      de=.5e5 * (ecp-ecm)
      dde= 1.e10 * ( (ecp-ec) + (ecm-ec) )
      gm=.25-.0682068*rs**2*(rs*dde-2*de)
      g0=.125/bi1(1.62903*sqrt(rs))**2
      a=.029
      b=.5625*gm-3.*(1.-g0)/64. -16.*a/15.
      c=-.75*gm+.5625*(1.-g0)-3.2*a
      b1=b+8*a/3.
      fwvi=1./fwv
      rls=rs
      endif

      q=rk*fwvi
      q2=q**2
      g=c+q2*(b+a*q2)+(-c+q2*(b1+q2*a))*(ff-1.)
      return
      end
