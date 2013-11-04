      subroutine checkin(n,a,b,it)
      include 'mach.p'
      integer n,it,i,j,k
      real*8 a,b,cmax,c
      dimension a(n,n),b(n,n)
      cmax=0.d0
      do 10 i=1,n
      do 10 j=1,n
      c=0.d0
      if(i.eq.j) c=-1.d0
      do 11 k=1,n
11    c=c+a(k,i)*b(k,j)
      cmax=max(cmax,abs(c))
      if(abs(c).gt.1.e-3) then
        write (*,*) ' checkin error ',it, i , j, c
!       stop
      endif
10    continue
      return
      end
