
      function ismax(n,x,incx)
      include 'mach.p'
      integer ismax,n,incx,i,maxc
      real*8  x(1+(n-1)*incx),xmax
      xmax=x(1)
      maxc=0
        do 1 i=1,n-1
      if (x(1+i*incx).gt.xmax) then
        xmax=x(1+i*incx)
        maxc=i
      endif
1     continue
      ismax=maxc+1
      return
      end
