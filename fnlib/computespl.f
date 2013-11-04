       subroutine computespl(rpos,w,m,mp,maxn,nknots,s,r,t,delta,coul)
       implicit none
       integer i,alpha,ialpha,i1alpha,n,maxn,m,nknots,maxim,mp
       double precision r(0:nknots),s(0:mp,0:2*mp+1)
       double precision t(0:nknots*(m+1)-1)
       double precision rpos,delta,w
       double precision rpi,rpi1,rinv,hplus,hminus
       double precision im,im1,iwgtminus
       logical coul

!   // rpos: distance from the origin(INPUT)
!   // w: value of the function
!   // m: number of continuous derivative(INPUT)
!   // maxn:maximum power of r in the spline-basis
!   // nknots:number of knots in the spline grid(INPUT)
!   // s: coefficient matrix (Natoli-Ceperley)(INPUT)
!   // r: position of the knots(INPUT)
!   // t: coefficient weighting the splines(INPUT)
!   // delta: distance between knots(INPUT)
!   // coul: 1/r factor or not(INPUT)
 
       maxim=nknots*(m+1)-1
!   // NEIGHBOR KNOT ON THE LEFT
       i=int(rpos/delta)
       if (i.ge.nknots) then
        w=0.d0
        return
       endif
!   // NORMALIZED DISTANCES
       rpi=(rpos-r(i))/delta
       rpi1=(r(i+1)-rpos)/delta
!   // COULOMB INTERACTION
       rinv=1.d0
       if (coul) rinv=1.d0/rpos
!   // LAST INTERVAL
       iwgtminus=1
       if (i.eq.nknots-1) iwgtminus=0
!   // AUXILIARY INDEX
       im=i*(m+1)
       im1=mod((i+1)*(m+1),maxim) !this ensures array bounds are not exceeded

       w=0.d0
       do alpha=0,m
        hplus=0.d0
        hminus=0.d0
        do n=0,maxn
          hplus=hplus+s(alpha,n)*rpi**n
          hminus=hminus+s(alpha,n)*rpi1**n
        enddo
        hplus=delta**alpha*hplus*rinv
        hminus=(-delta)**alpha*hminus*rinv
        ialpha=im+alpha
        i1alpha=im1+alpha
        w=w+t(ialpha)*hplus+t(i1alpha)*hminus*iwgtminus
       enddo
       return
       end
