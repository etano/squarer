      subroutine samcg(rmid,drift,covar,drmax,nv,ndim,ifnow,rn,slog
     +   ,width,nsamp)
      include 'mach.p'
#define SIMPLE_SPRNG
#include "sprng_f.h"

      integer ndm,nv,ndim,ifnow,nvd,ifl,j,k,ia,kk,i,l,ks,nsamp,nslices
     .,ip

      parameter (ndm=33)
      real*8 drmax,width(nv),eps,ri,det,sum,slog0
      real*8 covar(nv,nv),drift(nv),rn(nv,nsamp),rmid(nv),chi(ndm)
     +     ,adr(ndm),temp(ndm),slog(nsamp)
       include 'cpbc.cm'
c samples new coordinates rn based on rmid, drift, and covariance
c slog is the log of the sampling probability (returned)
c if for ifnow=2 sample new position, ifnow.eq.1 backwards step-just find slog
c drmax is the maximum drift
c ndim is the spatial dimensionality,nv is ndim*nmovers
      if(nsamp.gt.1)stop
      if(nv.gt.ndm) then
         write (*,*)' ndm too small in samcg ',nv,ndm
         stop
      endif
      nvd=nv/ndim
      slog0=0.d0

c puts choleski decomposition of covar back into covar.
      ifl=0
      i=0
      do ip=1,nvd
      do l=1,ndim
       i=i+1
      do 4 j=1,i-1
      do 4 k=i,nv
4     covar(k,i)=covar(k,i)-covar(k,j)*covar(i,j)
       if(covar(i,i).le..1*width(ip).or.ifl.eq.1) then
        ifl=1
        covar(i,i)=width(ip)
        do 8 k=i+1,nv
8       covar(k,i)=0.d0
        else
c     if(covar(i,i).le.0.) then
c      write (*,*)'matrix in sampr not positive definite',i,ifnew,nv
c      write (*,*) ((covar(i1,i2),i1=1,nv),i2=1,nv)
c      stop
c     endif
      covar(i,i)=sqrt(covar(i,i))
      ri=1.d0/covar(i,i)
      do 7 k=i+1,nv
7     covar(k,i)=ri*covar(k,i)
      endif
      enddo
      enddo

c now compute determinant of covar for the sampling probability
      det=covar(1,1)
      do 210 i=2,nv
210   det=det*covar(i,i)
      slog0=slog0-log(det)
 
c limit drift to drmax
      do 310 ia=1,nv
310   adr(ia)=max(-drmax,min(drmax,drift(ia)))
 
      do ks=1,nsamp
      slog(ks)=slog0
      if(ifnow.eq.2) then
c this is a forward step: sample nv gaussian numbers
      do 53 l=1,nv
      temp(l)=0.d0
53    chi(l)=sqrt(-2.d0*log(sprng()))
      do l=1,nv
      chi(l)=chi(l)*cos(3.14159265358979d0*sprng())
      enddo
c finds temp where temp=covar*chi
      do 3 j=1,nv
      do 3 k=j,nv
3     temp(k)=temp(k)+covar(k,j)*chi(j)
      k=0
      do kk=1,nvd
      do l=1,ndim
      k=k+1
!     rn(k,ks)=fpbc(temp(k)+rmid(k)+adr(k),l)
      rn(k,ks)=temp(k)+rmid(k)+adr(k)
      if(rn(k,ks).gt. el2(l))rn(k,ks)=rn(k,ks)-ell(l)
      if(rn(k,ks).lt.-el2(l))rn(k,ks)=rn(k,ks)+ell(l)
      enddo
      enddo
      else
 
c this is for a reverse step
      k=0
      do kk=1,nvd
      do l=1,ndim
      k=k+1
!      temp(k)=fpbc(rn(k,1)-rmid(k)-adr(k),l)
      temp(k)=rn(k,1)-rmid(k)-adr(k)
      if(temp(k).gt. el2(l))temp(k)=temp(k)-ell(l)
      if(temp(k).lt.-el2(l))temp(k)=temp(k)+ell(l)
      enddo
      enddo
c finds chi where temp=covar*chi
      do 81 j=1,nv
      sum=temp(j)
      do 71 k=1,j-1
71    sum=sum-covar(j,k)*chi(k)
81    chi(j)=sum/covar(j,j)
      endif
 
c for either direction find the probability
c we neglect the images
      do 22 k=1,nv
22    slog(ks)=slog(ks)-.5d0*chi(k)*chi(k)
      enddo
      end
