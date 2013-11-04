      real*8 function delta_action(x,y,s)
! ehbs2m is hbar^2/(2 mu).  (relative mass) tau= imag time slice.
! 1d action for 2 particles relative distances at the two time slices of x and y, 
!  x and y are assumed to be multiplied by  z=1.d0/dsqrt(4.d0*ehbs2m*tau)  on input
!interaction g*delta(x)  g negative means attractive
!     s=g*dsqrt(tau/(4.d0*ehbs2m)) 
! the various cutoffs (below) are chosen to make the error smaller than 1.e-8
      implicit none
      real*8 x,y,s,r,a,a2,ef

      r=dabs(x)+dabs(y)
      if(r+s.gt.8.d0) then !expand erfc in powers of (r+s)^2
         if(x*y.gt.5.d0) then
               delta_action=0.d0
               return
         endif
         ef=1.d0
         if(x*y.gt.0.d0)  ef=ef*dexp(-4.d0*x*y)
         a=1.d0/(r+s)
         a2=0.5d0*a*a
         ef=ef*s*a*(1.d0+a2*(-1.d0+3.d0*a2*(1.d0-5.d0*a2)))! error is 6a^9<4e-8
         delta_action=-dlog(1.d0-ef)
      elseif(r+s.lt.-5.d0) then
        a=s*(r+r+s)+(x-y)**2
        if(a.gt.23.d0) then
           delta_action=-a-dlog(-3.5449077018d0*s)
        elseif(a.lt.-23.d0) then
           delta_action=0.d0
        else
           delta_action=-dlog(1.d0-s*dexp(a)*3.5449077018d0) !sqrt(4pi)
        endif
      else ! full expression
        delta_action=-dlog(1.d0
     &  -s*dexp(s*(r+r+s)+(x-y)**2)*derfc(r+s)*1.77245385091d0) !sqrt(pi)
      endif
      end
