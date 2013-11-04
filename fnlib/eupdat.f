      subroutine eupdat(n,a,abar,error)
      include 'mach.p'
c updates cumlative mean and error with present value, n is the block number
      integer n
      real*8 a,abar,error,asq,n1,ni
      if(n.le.0) then
        write (*,*) ' n too small in eupdat ',n
        stop
      elseif(n.eq.1) then
        error=0.d0
        abar=a
      else
      n1=dfloat(n-1)
      ni=1.d0/dfloat(n)
      asq=n1*(error**2*dfloat(n-2)+abar**2)
      abar=(n1*abar+a)*ni
      error=sqrt(abs((asq+a**2)*ni-abar**2)/n1)
      endif
      return
      end
