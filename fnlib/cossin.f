      subroutine cossin(x,sp,n2,nk,ndim,rkcomp)
      include 'mach.p'
c calculates nk sin and cos for a particle located at x
c and puts them in sp.
      integer n2,nk,ndim,i
      real*8 x(ndim),rkcomp(ndim,nk)
      real*4 dum,sp(n2,nk)

      if(ndim.eq.3) then
        do i=1,nk
          dum=rkcomp(1,i)*x(1)+rkcomp(2,i)*x(2)+rkcomp(3,i)*x(3)
          sp(1,i)=cos(dum)
          sp(2,i)=sin(dum)
        end do

      else if (ndim.eq.2) then
        do i=1,nk
          dum=rkcomp(1,i)*x(1)+rkcomp(2,i)*x(2)
          sp(1,i)=cos(dum)
          sp(2,i)=sin(dum)
        end do

      else if (ndim.eq.1) then
        do i=1,nk
          dum=rkcomp(1,i)*x(1)
          sp(1,i)=cos(dum)
          sp(2,i)=sin(dum)
        enddo
      endif

      end
