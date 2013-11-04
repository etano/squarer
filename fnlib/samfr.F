      subroutine samfr(rmid,nmovers,ndim,rt,width,slog,inow,nsamp)
      include 'mach.p'
#define SIMPLE_SPRNG
#include "sprng_f.h"

c sample a new position centered at rmid

      integer nmovers,ndim,inow,nslices,i,l,nsamp,ks
      real*8 rmid(ndim,nmovers),rt(ndim,nmovers,nsamp),width(nmovers)
     +,slog(nsamp),a1,a2,chi,w2,x1,x0,x2

      include 'cpbc.cm'

      if(nsamp.ne.1)stop !NOTE TRAP
      do ks=1,nsamp
      if(inow.eq.2) then
        do i=1,nmovers
         do l=1,ndim
c make sure compiler does these in right order
          a1=sprng()
          a2=sprng()
          chi=sqrt(-2.d0*width(i)*log(a1))*cos(3.14159265358979d0*a2)
          rt(l,i,ks)=fpbc(rmid(l,i)+chi,l)
         end do
        end do
      endif
      slog(ks)=0.d0

      do i=1,nmovers
       w2=.5d0/width(i)
       do l=1,ndim
c put in the two nearest images
        x1=fpbc(rt(l,i,1)-rmid(l,i),l)
        x0=x1*x1
        slog(ks)=slog(ks)-w2*x0
        x2=x1-sign(ell(l),x1)
        slog(ks)=slog(ks)+log(1.d0+exp(-w2*(x2*x2-x0)) )
       end do
      end do
      end do

      return
      end
