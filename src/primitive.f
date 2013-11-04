      subroutine primitive(nx,rv,pott,ntemp,taulow,ehbs2m,diag
     +,mx,ndim,md)
      implicit none
       integer i,itemp,nx,ntemp,mx,ndim,md
      real*8 rv(nx),pott(nx),diag(mx,md,ntemp),gfun,dr,dri,tau
     +,s,taulow,ehbs2m
c find laplacian of potential
      do i=1,nx
      gfun=0.
      if(i.ne.1.or.i.ne.nx) then
        dr=.5*(rv(i+1)-rv(i-1))
        dri=1./dr
        gfun=dri**2*(pott(i+1)+pott(i-1)-2*pott(i))
     +          +(ndim-1)*(pott(i+1)-pott(i-1))/(2*dr*rv(i))
      endif

      tau=taulow
        do itemp=1,ntemp
        s=tau*ehbs2m/6.
c note that we have zeroed this out (primitive approximation)
        s=0.d0

         diag(i,1,itemp)=tau*(pott(i)+.5*s*gfun)
         if(md.ge.2)diag(i,2,itemp)=s*gfun
         if(md.ge.3)diag(i,3,itemp)=(.5*s)*gfun
         tau=tau+tau
        enddo

      enddo
      return
      end
