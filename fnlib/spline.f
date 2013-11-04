      subroutine spline(f,n,m,ifderv)
      include 'mach.p'
c sets up spline table, function values f on the n gridpoints are input
c if ifderv=1 the first derivates are estimated; otherwise they are input
c note that these are the derivates times delta(x)

      real*8 f,tw,th
      integer n,m,i,ifderv
      dimension f(m,4)
      if(n.gt.m) then
        write (*,*) ' spline error ',n,m
        stop
      endif

      if(ifderv.eq.1) then
        tw=1.d0/12.d0
        th=1.d0/1.5d0
        do 1 i=3,n-2
          f(i,2)=tw*f(i-2,1)-th*f(i-1,1)+th*f(i+1,1)-tw*f(i+2,1)
1       continue
        f(2,2)=.5*(f(3,1)-f(1,1))
        f(n-1,2)=.5*(f(n,1)-f(n-2,1))
        f(1,2)=2.d0*f(2,2)-f(3,2)
ccc(mdj)        f(n,2)=2.d0*f(n-1,2)-f(n-2,2)
        f(n,2)=0.d0
      endif

c determine 2rd and 3th derivative by matching
      do 2 i=1,n-1
        f(i,3)=3.d0*(f(i+1,1)-f(i,1))-2.d0*f(i,2)-f(i+1,2)
        f(i,4)=f(i+1,2)+f(i,2)+2.d0*(f(i,1)-f(i+1,1))
2     continue
c we just use linear extrapolation-but these values will not be used
      f(n,3)=2.d0*f(n-1,3)-f(n-2,3)
      f(n,4)=2.d0*f(n-1,4)-f(n-2,4)
      return
      end
