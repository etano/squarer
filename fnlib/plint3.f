      subroutine plint3(r,n,c,i0,pn,np)
      implicit none
c integrates sin(r*x)*x**i for i=i0 to n-1+i0 and x from 0 to 1
c pn=solid angle /volume, multiply by (x-1)**np
c result goes into c
      integer n,i0,np,i
      real*8 c(n+np),r,pn,ri
      complex*8 ti,et,em

      if(abs(r).gt.1.d-10) then
         ri=1.d0/r
         ti=cmplx(0.d0,-ri)
         et=cmplx(sin(r)*ri,-cos(r)*ri)
         em=ti*(et-ti)
         do i=1,n+i0+np
            if(i.gt.i0)c(i-i0)=real(em)
            em=ti*(et-i*em)
         enddo
      else
          do i=1,n+np
            c(i)=1.d0/(i+1.d0+i0)
          enddo
      endif

      call pmult(c,n,np,pn)
      return
      end
