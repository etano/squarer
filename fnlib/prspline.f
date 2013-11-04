      subroutine prspline(f,n,m,up)
      include 'mach.p'
c prints out a check of the derivitives and functions

      real*8 f,up,csiv
      integer n,m,i

      dimension f(m,4)
c now scale the data and spline it
      csiv=up/(n-1.)
       do 20 i=1,n
20     f(i,2)=f(i,2)*csiv
c make sure that endpoint values are zero
       f(n,1)=0.
       f(n,2)=0.
c check derivative 
c     vcum=f(n,1) 
c     do 15 i=n-1,1,-1 
c     vcum=vcum-.5*f(i,2) 
c     if(f(i,1).ne.0.) then 
c        f(i,3)=vcum/f(i,1)-1. 
c      else  
c        f(i,3)=0. 
c      endif 
c15    vcum=vcum-.5*f(i,2) 
c ia=number to entries printed, ip=number of "derivatives" printed
c      ia=10 
c      ip=1 
c      is=n/ia 
c     write (6,*) ' relative error of integrated derivative'
c     write (6,50) (f(i*is,3),i=1,ia) 
50    format(5e14.6)
c     write (6,*) ' table entries ' 
c     write (6,50) ((f(i*is,k),i=1,ia),k=1,ip)
      call spline(f,n,m,2)
      return
      end
