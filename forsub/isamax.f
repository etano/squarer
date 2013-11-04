
      function isamax(n,x,incx)
      include 'mach.p'
      integer n,incx,isamax,i,maxc
      real*8 x(1+(n-1)*incx),xmax
      xmax=abs(x(1))
      maxc=0
        do 1 i=1,n-1
      if (abs(x(1+i*incx)).gt.xmax) then
        xmax=abs(x(1+i*incx))
        maxc=i
      endif
1     continue
      isamax=maxc+1
      return
      end
