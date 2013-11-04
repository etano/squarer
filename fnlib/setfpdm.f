      subroutine setfpdm(tau,nslices,hbs2m,fpdm,mfpdm,morder,ndim,ell)
      include 'mach.p'

      real*8 tau,hbs2m,fpdm,ell,pi,arg,dr,w,wi
      real*8 c0,c1,c2,c3,r,ex,ex1
      integer nslices,mfpdm,morder,ndim,k,l,i,n

      dimension fpdm(mfpdm,morder,ndim,nslices),ell(ndim)
      pi=3.1415926535d0
      if(hbs2m*tau*nslices.le.0.) then
         write (*,*)' bad argument in setfpdm ',hbs2m,tau,nslices
c        stop
         return
      endif

      do 10 k=1,nslices
       do 10 l=1,ndim
        arg=(.5d0*ell(l))**2/(hbs2m*tau*k)
        dr=ell(l)/(mfpdm-1)
  
        if(arg.gt.pi) then
c image expansion
          do 12 i=1,mfpdm
            fpdm(i,1,l,k)=0.d0
            fpdm(i,2,l,k)=0.d0
12        continue
          w=4.d0*hbs2m*k*tau
          wi=1.d0/w
          c0=1.d0/sqrt(pi*w)
          c1=-2.d0*wi*dr
          do 14 n=-4,4
           do 14 i=1,mfpdm
            r=(n-.5d0)*ell(l)+(i-1)*dr
            ex1=wi*r**2
            ex=c0*exp(max(-ex1,-4000.))
            fpdm(i,1,l,k)=fpdm(i,1,l,k)+ex
            fpdm(i,2,l,k)=fpdm(i,2,l,k)+c1*r*ex
14        continue

        else

c     eigenfunction expansion
          c1=1.d0/ell(l)
          do 24 i=1,mfpdm
            fpdm(i,1,l,k)=c1
            fpdm(i,2,l,k)=0.d0
24        continue
          do 20 n=1,4
            c0=pi*n*(mfpdm+1.d0)/(mfpdm-1.d0)
c bug fixed 11/17/92
            c1=2.*exp(max(-(n*pi)**2/arg,-4000.))/ell(l)
            c2=2*pi*n*dr/ell(l)
            c3=c1*c2
            do 20 i=1,mfpdm
            fpdm(i,1,l,k)=fpdm(i,1,l,k)+c1*cos(c2*i-c0)
            fpdm(i,2,l,k)=fpdm(i,2,l,k)-c3*sin(c2*i-c0)
20        continue

        endif

        call spline(fpdm(1,1,l,k),mfpdm,mfpdm,0)

10    continue
      return
      end
