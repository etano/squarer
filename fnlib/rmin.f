      real*8 function rmin(r,ndim,nf,grad)
      include 'mach.p'
c finds minimum distance between two particles
      integer ndim,nf,nslices,j,i,im,jm,l
      real*8 r(ndim,nf),s,sij,grad(ndim,nf),g
      include 'cpbc.cm'
      s=ndim*el2(1)**2

      if (ndim.eq.2) then
      do 20 i=1,nf
      do 20 j=i+1,nf
      sij=fabc(r(1,i)-r(1,j),1) +fabc(r(2,i)-r(2,j),2)
       if(sij.lt.s) then
        im=i
        jm=j
        s=sij
       endif
20     continue

      else ! three dimensions
      do 30 i=1,nf 
      do 30 j=i+1,nf 
      sij=fabc(r(1,i)-r(1,j),1) +fabc(r(2,i)-r(2,j),2)
     &   +fabc(r(3,i)-r(3,j),3) 
       if(sij.lt.s) then
        im=i
        jm=j
        s=sij
       endif
30    continue
      endif

      rmin=sqrt(.5d0*s) ! This is distance to a node
      s=1.d0/s
      do l=1,ndim
        do i=1,nf
          grad(l,i)=0.d0
        enddo
        g=(fpbc(r(l,im)-r(l,jm),l))*s ! this is logarithmic derivative
        grad(l,im)=g
        grad(l,jm)=-g
      enddo
      return
      end
