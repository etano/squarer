
      function sang(ndim)
      include 'mach.p'
c solid angle in ndim dimensions
      real*8 sang,pi
      integer ndim
      pi=3.1415 92653 58979
      sang=2*pi*(ndim-1)
      if(ndim.eq.1)sang=1.
      return
      end
