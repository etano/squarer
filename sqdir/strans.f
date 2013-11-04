      subroutine strans(gr,rv,nr,gk,rkv,nk,ndim)
      implicit none
      integer nr,nk,ndim,k,i
      real*8 gr(nr),gk(nk),rv(nr),rkv(nk),g0,wt,bessj0
c puts the  transform of gr into gk (from 0 to rv)
c     integral of cos(kr)*gr(r)     if ndim=1
c     integral of j0(kr) *gr(r)     if ndim=2
c     integral of sin(kr)*gr(r)     if ndim=3
 
      g0=0.d0
      if(ndim.eq.1)g0=.5d0*gr(1)*rv(1)
        do k=1,nk
           gk(k)= g0
        enddo
 
      do i=1,nr
       wt=.5d0*gr(i)
       if(i.eq.1) then
        wt=wt*rv(2)
       elseif(i.eq.nr) then
        wt=wt*(rv(nr)-rv(nr-1))
       else
        wt=wt*(rv(i+1)-rv(i-1))
       endif
 
       if(ndim.eq.1) then
          do  k=1,nk
             gk(k)=gk(k)+wt*cos(rv(i)*rkv(k))
          enddo
       elseif(ndim.eq.2) then
          do  k=1,nk
             gk(k)=gk(k)+wt*bessj0(rv(i)*rkv(k))
          enddo
       elseif(ndim.eq.3) then
          do  k=1,nk
             gk(k)=gk(k)+wt*sin(rv(i)*rkv(k))
          enddo
       endif
      enddo
      return
      end
