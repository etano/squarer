      subroutine mdlng(ndim,p,ell,e2,v)
c Computes the madelung constant for r**-p interaction in ndim dimensions.
      implicit none
      integer icount(3),ndim,i,l,limit
      real*8 ell(ndim),pi,eps,vol,alpha,g5,v1,v2,v,r2,rk2,p,p2
     .,e2,gk0,gk,gr0,gr,gamma
      pi=3.1415 92653 58979d0
c This is accuracy of the constant. (but check how accurate gammi is?)
      eps=1.d-9
      limit=1+sqrt(-log(eps)/pi)
      vol=1.d0  
      do l=1,ndim
          vol=vol*ell(l)
          icount(l)=-limit
      enddo

      alpha=pi*vol**(-2.d0/ndim)

      g5=gamma(.5d0*p)
      v1=-2*alpha**(.5d0*p)/(p*g5)
      v2=-2*pi**(.5d0*ndim)/((ndim-p)*vol*g5*alpha**(.5d0*(ndim-p)))
      v=v1+v2 

      do 2 i=1,(2*limit+1)**ndim
      r2=0.d0   
      rk2=0.d0
      do l=1,ndim
           r2=r2+(icount(l)*ell(l))**2
           rk2=rk2+(2*pi*icount(l)/ell(l))**2
      enddo
      if(r2.gt.0)  then
         p2=.5d0*p 
         call gammi(gr,p2,alpha*r2,gr0)
         call gammi(gk,.5d0*(ndim-p),rk2/(4.d0*alpha),gk0)
         v=v+pi**(.5d0*ndim)*(4.d0/rk2)**(.5d0*(ndim-p))*gk/(vol*gr0)
     +       +gr/(gr0*r2**(p2))
      endif

      do 5 l=1,ndim
      icount(l)=icount(l)+1
      if(icount(l).le.limit) go to 2
5     icount(l)=-limit

2     continue

c multiply by charge**2
      v=v*e2
      return
      end
