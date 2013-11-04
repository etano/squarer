      function dleak(ifl,y)
      implicit none
c 1-D probability of leaking across a node=pnx, dpnx=is tau derivative
c y=x1*x2/(lambda*tau) is the input. x1 and x2 are nodal distances are 2 time steps
c ifl=0 return action, ifl=1 return derivative
c exclude y large since there will be no contribution
c and y=0 since that means nodal distance=infinity
      real*8 dleak,y
      integer ifl

      if(y.gt.14.d0) then ! neglect contribution of 1e-6
         if(ifl.eq.0)dleak=1.d0
         if(ifl.eq.1)dleak=0.d0
      else
        if(y.le.0.d0) then
           write (*,*) 'problem in dleak ifl=',ifl,' arg=',y
!          stop
           if(ifl.eq.0) dleak=0.d0
           if(ifl.eq.1) dleak=1.d0
        endif
        if(ifl.eq.0)dleak=(1.d0-dexp(-y))
        if(ifl.eq.1)dleak=y/(dexp(y)-1.d0)
      endif
      end
