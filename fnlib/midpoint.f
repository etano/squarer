      subroutine midpoint(r1,r2,rmid,nmovers,ndim)
      include 'mach.p'
      integer nmovers,ipn,l,nslices,ndim
      real*8 del,rmid(ndim,nmovers),r1(ndim,nmovers),r2(ndim,nmovers)
      include 'cpbc.cm'
c for all movers construct mean positions 
      do ipn=1,nmovers
c  rmid=.5*(r1+r2)
        do l=1,ndim
          del=fpbc(r2(l,ipn)-r1(l,ipn),l)
          rmid(l,ipn)=fpbc(r1(l,ipn)+.5d0*del,l)
        enddo
      enddo
      end
