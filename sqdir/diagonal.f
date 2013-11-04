      subroutine diagonal(nx,ndim,tau,ehbs2m,rv,lmax,uold,lpos,nl,mx
     + ,diag,ml,mx2,pott,md)
c this determines the diagonal density matrix by summing over partial waves.

      implicit none
      integer ldim,ndim,nx,lmax,nl,lul,i,l,mx,isym,id,k
     +,ml,mx2,md
      parameter (ldim=128)
      real*8 ehbs2m,tau,cp,coeff,pi4,etau,cutl,s1,s2,s3,ss(3),coefr,c
     +,rhol,rhol0,drhol0,eul,gs,corr
      real*8  rv(nx),cl(0:ldim),diag(mx,3),f(0:ldim,3),free(0:ldim,2)
     +,uold(0:ml,mx2,3),pott(nx)
      integer lpos(nx)

c initialize partial wave components
      call pwc(ndim,ldim,cl)
      etau=ehbs2m*tau
      cp=.5d0*ndim/tau
      pi4=12.566 37061 43591d0
      coeff=sqrt(pi4*etau)**ndim
      cutl=-log(1.d-10)
      lmax=1

      do i=1,nx
 
        coefr=coeff/rv(i)**(ndim-1)
        lul=min0(nl,2+int(sqrt(cutl*rv(i)**2/etau)) ) 
        call fdenmt(rv(i),rv(i),free,etau,ldim,ehbs2m,lul,ndim)
 
        s1=0.d0
        s2=1.d0
        s3=0.d0
        do k=1,md
        ss(k)=0.d0
        enddo

        id=isym(i,i)
 
        do l=0,lul
           eul=exp(-uold(l,id,1))
 
           if(s2.gt.1.d-5) then
            corr=eul*s2
            f(l,1)=s1+corr
            do k=2,md
             f(l,k)=ss(k)-uold(l,id,k)*corr-s3*eul
            enddo
           else
c we cut off the corrections at lpos
            lpos(i)=l-1
            go to 7229
           endif
  
           c=cl(l)*coefr
           rhol0=c*free(l,1)
           rhol=rhol0*eul
           s1=s1+rhol
           s2=s2-rhol0

           drhol0=rhol0*(free(l,2)+cp)
           s3=s3+drhol0
           do k=2,md
            ss(k)=ss(k)+drhol0*eul-uold(l,id,k)*rhol
           enddo
        enddo
 
c lpos is the last value of l giving reasonable values for sum rule
        lpos(i)=lul

7229    gs=f(lpos(i),1)

        diag(i,1)=-log(gs)
        do id=2,md
        diag(i,id)=-f(lpos(i),id)/gs
        diag(i,id)=-f(lpos(i),id)/gs
        enddo
 
c convergence criterion for density matrix
        do l=0,lpos(i)-1
           if(abs(f(l,1)-gs).gt.1.e-4*gs.and.diag(i,1).lt.7.)
     &        lmax=max0(lmax,l)
         enddo

      enddo

      return
      end
