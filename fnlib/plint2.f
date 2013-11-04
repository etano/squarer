      subroutine plint2(r,n,c,i0,pn,np)
      implicit none
      integer n,i0,np,mbf,ier,nbf,nn,i,k
      parameter (mbf=6000)
      real*8 c(20),bess(mbf+2),r,pn,fact,con,rat,term
c integrates besj0(r*x)*x**i for i from i0-1 to n+i0-2
c using formula 11.1.1 from Abramowitz and Stegun

      if(r.ge.1.d-5) then
c determine mbf bessel functions
       nbf=max(100,6*int(r))
       if(nbf.gt.mbf) then 
         write (*,*) ' danger not enough space in plint2'
         nbf=mbf
       endif
       call mmbsjn(r,nbf,bess,ier)
!      if(ier.ne.0) write (*,*) 'problem in plint2 ',ier
       fact=1.d-10
       con=2.d0/(exp(1.d0)*r)
      endif

      do 1 i=1,n+np
      nn=i+i0-2

      if(r.lt.1.d-5) then
      c(i)=1.d0
      else

      c(i)=bess(2)
      rat=1.d0
      do  k=4,nbf,2
      rat=rat*dfloat(-nn+k-4)/dfloat(nn+k)
      term=rat*bess(k)*dfloat(k-1)
      c(i)=c(i)+term
c stopping criterea: next term will be small relative to total.
      if(abs(rat).lt.(con*k)**k*fact) go to 5
      enddo
c     write (*,*) ' loop does not terminate in plint ',r,i,nbf
5     c(i)=2.d0*c(i)/r
      endif

1      c(i)=c(i)/dfloat(i+i0)
      call pmult(c,n,np,pn)
!     if(r.gt.81.5.and.r.lt.82.4) then
!       write (55,'(5e14.5)')(bess(i),i=1,n)
!     endif
      return
      end
