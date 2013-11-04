       subroutine splint3D(maxn,r,delta,k,vol,coul,ddplus,ddminus)
       implicit none
       integer maxn,i,n
       double precision r,ri,delta,k,vol
       double complex expri(-1:1),ee(0:maxn+1,-1:1),imunit
       double precision ddplus(0:maxn),ddminus(0:maxn)
       double precision auxplus(0:maxn+2),auxminus(0:maxn+2)
       double precision dnorm,inv1,inv2,inv3,pi
       logical coul
      parameter  (pi=3.14159265359d0)
       imunit=dcmplx(0.d0,1.d0)

!   // k=0 MUST BE TREATED DIFFERENTLY
       if(k.gt.1.d-8) then

!   // COMPUTE exp(ikr_i),exp(ikr_{i+1}),exp(ikr_{i-1})
       do i=-1,1
        ri=r+i*delta
        expri(i)=exp(imunit*k*ri)
       enddo

!   // RECURSIVELY COMPUTE E_{ikn}^{+-}
!   // E=ee; i=+,- 
       do i=-1,1,2
        ee(0,i)=-imunit/k*(expri(i)-expri(0))
        do n=1,maxn+1
         ee(n,i)=-imunit/k*(i**n*expri(i)-dble(n)/delta*ee(n-1,i))
        enddo
       enddo

!   // DEFINE THE VALUES OF D_{ikn}^{+-} COMBINING E_{ikn}^{+-}
!   // D_{ikn}^{+}=ddplus;D_{ikn}^{-}=ddminus
       dnorm=4.d0*pi/(k*vol)
       if (.not.coul) then
        do n=0,maxn
         ddplus(n)=dnorm*(delta*dimag(ee(n+1,1))+r*dimag(ee(n,1)))
         ddminus(n)=-dnorm*(delta*dimag(ee(n+1,-1))+r*dimag(ee(n,-1)))
        enddo
       else
        do n=0,maxn
         ddplus(n)=dnorm*dimag(ee(n,1))
         ddminus(n)=-dnorm*dimag(ee(n,-1))
        enddo
       endif

       else ! THE CASE OF k=0

       dnorm=4.d0*pi/vol*delta
       
       if (.not.coul) then
        do n=0,maxn+2
         auxplus(n)=1.d0/dble(n+1)
         auxminus(n)=auxplus(n)*(-1.d0)**(n+1)
        enddo
        do n=0,maxn
         ddplus(n)=dnorm*(delta**2*auxplus(n+2)+r**2*auxplus(n)
     &   +2*delta*r*auxplus(n+1))
         ddminus(n)=-dnorm*(delta**2*auxminus(n+2)+r**2*auxminus(n)
     &   +2*delta*r*auxminus(n+1))
        enddo
       else
        do n=0,maxn
         inv1=1.d0/dble(n+1)
         inv2=1.d0/dble(n+2)
         ddplus(n)=dnorm*(r*inv1+delta*inv2)
         ddminus(n)=dnorm*(r*inv1-delta*inv2)*(-1.d0)**(n+2)
        enddo
       endif

       endif

       return
       end

