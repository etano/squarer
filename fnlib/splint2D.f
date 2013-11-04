       subroutine splint2D(maxn,r,delta,k,vol,coul,ddplus,ddminus)
       implicit none
       integer maxn,i,n,l
       real*8 r,ri,delta,k,vol
       real*8 J0,J1,sumplus,summinus,dnorm,x
       real*8 ddplus(0:maxn),ddminus(0:maxn)
       real*8 bincoeff(0:maxn,0:maxn),int_j0rl(0:maxn+1,-1:1)
       real*8 h0,h1,pi,bessj0,bessj1
       logical coul
      parameter  (pi=3.14159265359d0)

!   // k=0 MUST BE TREATED DIFFERENTLY
       if(k.gt.1.d-8) then

!   // COMPUTE binominal coeff bincoeff(n,l)=(n l) and int_j0rl(l)=int r^l bessjo(r)
       do n=0,maxn
          bincoeff(n,0)=1.d0
          do l=1,n
             bincoeff(n,l)=(n-l+1.d0)/(1.d0*l)*bincoeff(n,l-1)
          end do
       end do
       do i=-1,1
        ri=r+i*delta
        x=k*ri
        J0=bessj0(x)
        J1=bessj1(x)
!
!            int dx     J0(x) = x J0(x) + pi/2 x ( H0(x) J1(x)- H1(x) J0(x) )
!            int dx x   J0(x) = x J1(x)
!
!            int dx x^l J0(x) = x^l J1(x) +(l-1) x^(l-1) J0(x) - (l-1)^2 int dx x^(l-2) J0(x)
!
!            int_j0rl(l,.)=k^(-l) int dx x^l J0(x)    (x=k*r_i)
!
        int_j0rl(0,i)=(x*J0+pi/2.d0*x*(h0(x)*J1-h1(x)*J0))
        int_j0rl(1,i)=ri*J1
        do l=2,maxn+1
           int_j0rl(l,i)= ri**l*J1+(l-1.d0)/(1.d0*k)*ri**(l-1)*J0
     .                    -((l-1.d0)/(1.d0*k))**2*int_j0rl(l-2,i)
        end do
       enddo


!   // CALCULATE THE VALUES OF D_{ikn}^{+-} USING the binominal expression and the interals
!   // D_{ikn}^{+}=ddplus;D_{ikn}^{-}=ddminus

       dnorm=2.d0*pi/(k*vol)
       if (.not.coul) then
        do n=0,maxn
         sumplus=0.d0
         summinus=0.d0
         do l=0,n
            sumplus=sumplus+bincoeff(n,l)*(-r)**l
     .                     *(int_j0rl(n+1-l,1)-int_j0rl(n+1-l,0))
            summinus=summinus+bincoeff(n,l)*(-r)**l
     .                     *(int_j0rl(n+1-l,-1)-int_j0rl(n+1-l,0))
         end do
         ddplus(n)=dnorm/delta**n*sumplus
         ddminus(n)=-dnorm/delta**n*summinus
        enddo
       else
        do n=0,maxn
         sumplus=0.d0
         summinus=0.d0
         do l=0,n
            sumplus=sumplus+bincoeff(n,l)*(-r)**l
     .                     *(int_j0rl(n-l,1)-int_j0rl(n-l,0))
            summinus=summinus+bincoeff(n,l)*(-r)**l
     .                     *(int_j0rl(n-l,-1)-int_j0rl(n-l,0))
         end do
         ddplus(n)=dnorm/delta**n*sumplus
         ddminus(n)=-dnorm/delta**n*summinus
        enddo
       endif

       else ! THE CASE OF k=0

       
       if (.not.coul) then
        dnorm=2.d0*pi/vol*delta**2
        do n=0,maxn
         ddplus(n)=dnorm*(1.d0/(n+2.d0)+r/(delta*(n+1.d0)))
         ddminus(n)=dnorm*(-1)**(n+1)*(1.d0/(n+2.d0)-r/(delta*(n+1.d0)))
        enddo
       else
        dnorm=2.d0*pi/vol*delta
        do n=0,maxn
         ddplus(n)=dnorm/(n+1.d0)
         ddminus(n)=dnorm*(-1)**n/(n+1.d0)
        enddo
       endif

       endif

       return
       end

