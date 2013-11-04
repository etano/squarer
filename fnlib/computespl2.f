       subroutine computespl2(rpos,w,dw,ddw,d3w,m,mp
     .          ,maxn,nknots,s,r,t,delta)
       implicit none
       logical coul
       integer i,alpha,ialpha,i1alpha,n,maxn,m,nknots,maxim,mp
       double precision r(0:nknots),s(0:mp,0:2*mp+1)
       double precision t(0:nknots*(m+1)-1)
       double precision rpos,delta,w
       double precision rpi,rpi1,rinv,hplus,hminus
       double precision im,im1,iwgtminus
       double precision dw,ddw,dhplus,dhminus,ddhplus,ddhminus
       double precision d3w,d3hplus,d3hminus

!   // rpos: distance from the origin(INPUT)
!   // w: value of the function
!   // dw: first derivative
!   // ddw: second derivative 
!   // m: number of continuous derivative(INPUT)
!   // maxn:maximum power of r in the spline-basis
!   // nknots:number of knots in the spline grid(INPUT)
!   // s: coefficient matrix (Natoli-Ceperley)(INPUT)
!   // r: position of the knots(INPUT)
!   // t: coefficient weighting the splines(INPUT)
!   // delta: distance between knots(INPUT)
! if coul=true divide values by r, otherwise do not
 
       maxim=nknots*(m+1)-1
!   // NEIGHBOR KNOT ON THE LEFT
       i=int(rpos/delta)
       if (i.ge.nknots) then
        w=0.d0
        dw=0.d0
        ddw=0.d0
        d3w=0.d0
        return
       endif
!   // NORMALIZED DISTANCES
       rpi=(rpos-r(i))/delta
       rpi1=(r(i+1)-rpos)/delta
!   // LAST INTERVAL
       iwgtminus=1
       if (i.eq.nknots-1) iwgtminus=0
!   // AUXILIARY INDEX
       im=i*(m+1)
       im1=mod((i+1)*(m+1),maxim) !this ensures array bounds are not exceeded

       w=0.d0
       dw=0.d0
       ddw=0.d0
       d3w=0.d0
       do alpha=0,m
        hplus=0.d0
        hminus=0.d0
        dhplus=0.d0
        dhminus=0.d0
        ddhplus=0.d0
        ddhminus=0.d0
        d3hplus=0.d0
        d3hminus=0.d0
        do n=0,maxn
          hplus=hplus+s(alpha,n)*rpi**n
          hminus=hminus+s(alpha,n)*rpi1**n
          if(n.gt.0) then
            dhplus=dhplus+s(alpha,n)*n*rpi**(n-1)
            dhminus=dhminus+s(alpha,n)*n*rpi1**(n-1)
            if(n.gt.1) then
              ddhplus=ddhplus+s(alpha,n)*n*(n-1)*rpi**(n-2)
              ddhminus=ddhminus+s(alpha,n)*n*(n-1)*rpi1**(n-2)
              if(n.gt.2) then
                d3hplus=d3hplus+
     .             s(alpha,n)*n*(n-1)*(n-2)*rpi**(n-3)
                d3hminus=d3hminus
     .              +s(alpha,n)*n*(n-1)*(n-2)*rpi1**(n-3)
              end if
            end if
          end if
        enddo
        hplus=delta**alpha*hplus
        hminus=(-delta)**alpha*hminus
        dhplus=delta**(alpha-1)*dhplus
        dhminus=(-delta)**(alpha-1)*dhminus
        ddhplus=delta**(alpha-2)*ddhplus
        ddhminus=(-delta)**(alpha-2)*ddhminus
        d3hplus=delta**(alpha-3)*d3hplus
        d3hminus=(-delta)**(alpha-3)*d3hminus

        ialpha=im+alpha
        i1alpha=im1+alpha
        w=w+t(ialpha)*hplus+t(i1alpha)*hminus*iwgtminus
        dw=dw+t(ialpha)*dhplus+t(i1alpha)*dhminus*iwgtminus
        ddw=ddw+t(ialpha)*ddhplus+t(i1alpha)*ddhminus*iwgtminus
        d3w=d3w+t(ialpha)*d3hplus+t(i1alpha)*d3hminus*iwgtminus
       enddo
!      if(coul) then
!        w=w/rpos
!        dw=(dw-w)/rpos
!        ddw=(ddw-2*dw)/rpos
!        d3w=(d3w-3*ddw)/rpos
!      endif
       end
