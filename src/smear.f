      subroutine smear(ndim,ntemp,ehbs2m,taulow,diag,rv,nx,mx,md)
c now determine sampling tables by convoluting with a Gaussian
c 1d potentials are not quite right
      implicit none
      integer nkv,ndim,ntemp,mtemp,nx,mx,mxx,i,l,k,id,j,nneg,md
c nkv are total number of k-values in fourier transforms for sampling tables
      parameter (nkv=2000,mxx=1000,mtemp=10)
      real*8 gk(nkv),rkv(nkv),diag(mx,md,ntemp),gfun(mxx)
     +,samtab(mxx,2,mtemp,2),f1,f2,a1,a2,sang
     +,ftnorm,pi,wdt,ehbs2m,taulow,rv(mx),u,dk,gfun0,wtt,delta,gwt

c danger for 1-d functions that do not go to 0 at large x or are not symmetric.
      if(ehbs2m.le.0.) return
      if(taulow.le.0.)stop
      if(ndim.lt.1.or.ndim.gt.3)stop
      if(nx.le.1)stop
      if(ntemp.le.0)stop
      if(nx.gt.mxx) then
        write (*,*)' mxx in smear too small '
        stop
      endif
      if(ntemp.gt.mtemp) then
        write (*,*)' mtemp in smear too small'
        stop
      endif
      pi=3.14159265358979d0
      if(ndim.eq.1)ftnorm=1.d0/(2.d0*pi)
      if(ndim.eq.2)ftnorm=1.d0
      if(ndim.eq.3)ftnorm=2.d0/pi
      write (*,*)' beginning smear'

      do 3000 l=1,ntemp
      wdt=.5d0*ehbs2m*taulow*2**(l-1)

c first make a table of gfun(r)=exp(-u(r))-1
       do i=1,nx
        u=diag(i,1,l)
        if(u.lt.300.) then
          gfun(i)=exp(-u)-1.d0
        else
         gfun(i)=-1.d0
        endif
        if(ndim.ge.2) gfun(i)=gfun(i)*rv(i)
       enddo

c make a uniform grid in k-space
c gaussian will drop by a factor of 1.e-9 at edge of grid
       dk=sqrt(-2.d0*log(1.d-9)/wdt)/nkv
       do  k=1,nkv
         rkv(k)=k*dk
       enddo

c take fourier transform into gk
       call strans(gfun,rv,nx,gk,rkv,nkv,ndim)

       if(ndim.eq.2) then
        do  k=1,nkv
         gk(k)=rkv(k)*gk(k)
        enddo
       endif

c j=1 is for mover-nonmover , j=2 is for mover-mover
       do 3070 j=1,2
c multiply by a gaussian to convolve
        do k=1,nkv
           gk(k)=gk(k)*exp(-.5d0*wdt*rkv(k)**2)
        enddo
c transform back
        call strans(gk,rkv,nkv,gfun,rv,nx,ndim)
c normalize 
        if(ndim.eq.3) then
         do i=1,nx
          gfun(i)=ftnorm*gfun(i)/rv(i)
         enddo
        elseif(ndim.le.2) then
         do i=1,nx
          gfun(i)=ftnorm*gfun(i)
         enddo
        endif
c diagnostic-gfun will be smeared version of exp(-u)-1 
c       do i=1,nx
c         write (89,89) rv(i),gfun(i),j,l
c89       format(2f12.5,2i5)
c       enddo

c for small r find value at r=0 and interpolate from that
        wtt=.5d0/(j*wdt)
        gwt=.5d0*sang(ndim)*(wtt/pi)**(.5d0*ndim)
        gfun0=0.d0
        do i=1,nx-1
          if(i.eq.1) then
            delta=gwt*rv(2)
          else
            delta=gwt*(rv(i+1)-rv(i-1))
          endif
        gfun0=gfun0+delta*rv(i)**(ndim-1)*exp(-diag(i,1,l)-wtt*rv(i)**2)
        enddo
        gfun0=-log(max(gfun0,1.d-12))

        do i=nx,1,-1
         if(gfun(i).le.-.999) go to 707
        enddo
        i=-1
707     nneg=i
        if(nneg.ge.nx-2)stop
        f1=(-log(gfun(nneg+2)+1.d0)-gfun0)/rv(nneg+2)**2
        f2=(-log(gfun(nneg+3)+1.d0)-gfun0)/rv(nneg+3)**2
        a2=(f2-f1)/(rv(nneg+3)**2-rv(nneg+2)**2)
        a1=f1-a2*rv(nneg+2)**2

        do i=1,nx
         if(i.le.nneg+1)then
c now interpolate assuming gfun=gfun0+ cnst*r**2
          gfun(i)=gfun0+rv(i)**2*(a1+a2*rv(i)**2)
         else
          gfun(i)=-log(1.+gfun(i))
         endif
        enddo

        if(ndim.eq.1) then
         do i=2,nx-1
c this is first derivative for the drift -wdt*(dw/dr)
          samtab(i,j,l,1)=-wdt*(gfun(i+1)-gfun(i-1))/(rv(i+1)-rv(i-1))
c this is second derivative for the covariance matrix
          samtab(i,j,l,2)=-wdt**2*( gfun(i+1)+gfun(i-1)-2*gfun(i) )
     +                     /(.5d0*(rv(i+1)-rv(i-1)) )**2
         enddo

        else
         do i=2,nx-1
c this is first derivative for the drift -wdt*(1/r)*(dw/dr)
          samtab(i,j,l,1)=-wdt*(gfun(i+1)-gfun(i-1))
     +                     /((rv(i+1)-rv(i-1))*rv(i))
c this is second derivative for the covariance matrix
c    -wdt**2 (1/r)*d/dr[ (1/r)*(dw/dr) ]
          samtab(i,j,l,2)=wdt*(-samtab(i,j,l,1)-wdt*(gfun(i+1)
     &                     +gfun(i-1)-2*gfun(i) )
     &                    /(.5d0*(rv(i+1)-rv(i-1)))**2  )  /rv(i)**2
         enddo
        endif
c end points
        do id=1,2
         samtab(1,j,l,id)=samtab(2,j,l,id)
         samtab(nx,j,l,id)=0.d0
        enddo
3070   continue

3000  continue

      do id=1,2
       write (1,*)'RANK 3 ',nx,' 2 ',ntemp
       write (1,*)'LABEL 1 r'
       write (1,*)'LABEL 2 movers/nonmovers'
       write (1,*)'LABEL 3 time'
       write (1,889)  3 ,taulow,taulow*2**ntemp
       write (1,*)'BEGIN sampling table  ',md+id
        do 3510 l=1,ntemp
        do 3510 j=1,2
3510    write (1,888)(samtab(i,j,l,id),i=1,nx)
889     format(' GRID ',i3,' LOG ',2e17.9)
888     format(5e17.9)
      enddo
      return
      end
