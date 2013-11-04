
      function isamin(n,x,incx)
      include 'mach.p'
      integer n,incx,isamin,i,minc
      real*8  x(1+(n-1)*incx),xmin
      xmin=abs(x(1))
      minc=0
        do 1 i=1,n-1
      if (abs(x(1+i*incx)).lt.xmin) then
        xmin=abs(x(1+i*incx))
        minc=i
      endif
1     continue
      isamin=minc+1
      return
      end
