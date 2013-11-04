      subroutine fitpn(n,v,rk,wt,nk,nf,b0,b,hs,i0,vol,a,ndim,ifcon,np
     +   ,vmad)
      implicit none
c version of 4/23/93 DMC
      integer n,nk,nf,i0,ndim,ifcon,np,mp,j,info,m,i2,i,k
      real*8 v(nk),rk(nk),wt(nk),b(n),hs(n,n+2),vmad,b0,vol,a,pi,pn
     .,beta1,beta2,anp,vmad1,vmad2,chisq,vv,vkold
c fits best spherical polynomial b to v(k); 
c first term is r**(i0-1). b0 is the value of r**(-1+2*i0) term.
c a is the cutoff in r. 
c ifcon=1 use the constraint b0 otherwise not
c Results are returned: in b, v, vmad
      pi=3.1415 92653 58979d0
      pn=dfloat(2*(ndim-1))*pi*a**ndim/vol
c np is the number of factors of (r-a) we multiply by
c m is the total number of free parameters
      m=n-np
      if(ifcon.ne.1) then
       mp=m
       i2=1
       beta1=0.d0
       beta2=0.d0
      else
c one more equation because of cusp condition
       mp=m-1
       i2=2
c We are setting up cusp condition constraints
       anp=dfloat(1-2*mod(np,2))
       if(i0.ne.0) then
        beta1=-anp*b0*a/float(np)
        beta2=1.d0/float(np)
       else  
         beta1=anp*b0/a
         beta2=0.d0
       endif 
      endif

c     write (6,*) ' beginning fitpn n nk nf i0= ',n,nk,nf,i0
c     write (6,*) ' ifcon ndim  ',ifcon,ndim
c     write (6,*) '  vol b0  ',vol,b0
c     write (6,*) ' a pn np m mp i2 beta ',a,pn,np,m,mp,i2,beta1,beta2
c zero right and left hand side of fitting equation
      do 1 i=1,n+2
      do 1 j=1,n
1     hs(j,i)=0.d0
      chisq=0.d0
c go over k values larger than those explicitly used
      do 20 k=nf+2,nk
c the matrix elements are different in 2 and 3 dimensions
      if(ndim.eq.3)call plint3(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      if(ndim.eq.2)call plint2(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      vv=v(k)-beta1*hs(1,m+2)
      chisq=chisq+wt(k)*vv**2
c for the derivative constraint
      hs(2,m+2)=hs(2,m+2)+beta2*hs(1,m+2)
c add to right hand side
      do 20 i=i2,m
              hs(i,m+1)=hs(i,m+1)+wt(k)*vv*hs(i,m+2)
c add to left hand side
              do j=i2,m
                    hs(j,i)=hs(j,i)+wt(k)*hs(i,m+2)*hs(j,m+2)
              enddo
20    continue

c invert right hand side
!     call spoco(hs(i2,i2),n,mp,rcond,hs(1,m+2),info)
      call spofa(hs(i2,i2),n,mp,info)
      if(info.ne.0) then
         write (*,*) ' trouble in fitpn info=',info
         stop
      endif
c     write (*,*) n,' condition number ',rcond

c make a spare copy of right hand side
      do i=i2,m
         b(i)=hs(i,m+1)
      enddo

c solve linear equations
      call sposl(hs(i2,i2),n,mp,hs(i2,m+1))

      do i=i2,m
         chisq=chisq-b(i)*hs(i,m+1)
      enddo
      if(chisq.gt.0) then
          write (6,*) 'rms fitting error=',sqrt(chisq)
       else
          write (6,*) ' chisq negative ? ',chisq
       endif

c this is cusp constraint
      if(ifcon.eq.1)hs(1,m+1)=beta1+beta2*hs(2,m+1)

c subtract effect of short range potential on fourier components
      do k=1,nk
       vkold=v(k)
           if(ndim.eq.3)call plint3(rk(k)*a,m,hs,i0,pn,np)
           if(ndim.eq.2)call plint2(rk(k)*a,m,hs,i0,pn,np)
           do i=1,m
              v(k)=v(k)-hs(i,m+1)*hs(i,1)
        enddo
        if(k.ge.212.and.k.le.213) then
        write (55,*) k,rk(k)*a,v(k),vkold
        write (55,'(5e12.5)') (hs(i,1),i=1,m)
        endif
      enddo

      do i=1,m
          b(i)=hs(i,m+1)
      enddo

      write(6,7) (b(i),i=1,m)
7     format('s.r. polynomial=',4e14.6)

c this is the derivative of the polynomial at x=0 for the self-interaction.
      vmad1=(-1.d0)**np*(b(2)-np*b(1))
      vmad2=0.d0
      do k=1,nk
         vmad2=vmad2+wt(k)*v(k)
!         write (1,*) k,wt(k),v(k),vmad2
      enddo
c the madelung constant is just for testing
      write(*,*)' computed madelung constant=',vmad1+vmad2,vmad1,vmad2
      write(*,*)'      gaussian charge value=',vmad, 
     + '  difference ', vmad-vmad2-vmad1
c return just the derivative at the origin.
      vmad=vmad-vmad2
!     vmad=vmad1 ! alternative because vmad2 has error?

      return
      end
