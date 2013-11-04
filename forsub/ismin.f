      function ismin(n,x,incx)
      include 'mach.p'
      integer ismin,n,incx,i,minc
      real*8 x(1+(n-1)*incx),xmin
      xmin=x(1)
      minc=0
        do 1 i=1,n-1
      if (x(1+i*incx).lt.xmin) then
        xmin=x(1+i*incx)
        minc=i
      endif
1     continue
      ismin=minc+1
      return
      end
