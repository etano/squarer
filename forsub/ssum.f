      function ssum(n,sx,incx)
      include 'mach.p'
      real*8 ssum,sx
      integer n,incx,i
      dimension sx(1+(n-1)*incx)
      ssum=0.
      do 1 i=0,n-1
      ssum=ssum+sx(1+i*incx)
1     continue
      return
      end
