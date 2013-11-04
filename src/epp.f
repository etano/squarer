      subroutine epp(q,rs,r,phir)
      implicit real*8 (a-h,o-z)
c finds the interaction between 2 protons in an electron gas
c rs is the electon gas density, the potential phir will be computed
c on the grid r. r is in atomic units. k in units of rs.
      parameter (nk=100000,dk=.01)
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
      end
