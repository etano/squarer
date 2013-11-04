      subroutine fdenmt(ri,rj,d,et,ldim,ehbs2m,lu,ndim)
      implicit none
c modify to calculate the derivative of the logarithm wrt to et
      integer ldim,lu,ndim,l,nout
      real*8 ri,rj,et,ehbs2m,pi,eps,wi,arg,eti,cc,a
      real*8 d(0:ldim,2)
      save pi,eps
      data pi/3.14159265358979d0/
ccc      data eps/1.e-7/
ccc adjust eps for a few more digits -- hopefully, see vmsbi
c
      data eps/1.d-10/
ccc
      wi=.25d0/et
      arg=2.d0*wi*ri*rj
      eti=ehbs2m/et

      if(ndim.eq.1)then 

       arg=wi*(ri-rj)**2
       d(0,1)=exp(-arg)/sqrt(4.d0*pi*et)
       d(0,2)=eti*(-.5d0+arg)

      elseif(ndim.eq.2) then
c slatec works better!! but now it has underflow problems
       call besi(arg,0.d0,2,lu+1,d,nout)
       if(nout.gt.0) then
!      write (*,*)' problem in besi',nout,lu,arg,d(0,1)
c numerical recipes desn't work? 1/22/93
       call bessi(lu,arg,d,eps)
       endif
       cc=sqrt(ri*rj)*exp(-(ri-rj)**2*wi)/(2.d0*et)
       do l=0,lu
        d(l,1)=d(l,1)*cc
       enddo
       a=-1.d0+(ri**2+rj**2)*wi
       d(0,2)=eti*(a-arg*d(1,1)/d(0,1))
       do l=1,lu
        d(l,2)=eti*(a+l-arg*d(l-1,1)/d(l,1))
       enddo

      elseif(ndim.eq.3) then

       d(0,1)=(exp(-(ri-rj)**2*wi)-exp(-(ri+rj)**2*wi))/sqrt(4.d0*pi*et)
       a=-.5d0+(ri**2+rj**2)*wi
       d(0,2)=eti*(a-arg/tanh(arg))
       if(lu.gt.0) then
        call vmsbi(arg,d(0,1),eps,lu+1)
        do l=1,lu
         d(l,2)=eti*(a+l-arg*d(l-1,1)/d(l,1))
        enddo
       endif

      endif

      return
      end
