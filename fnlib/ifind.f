      function ifind(a,b,n)
      include 'mach.p'
c finds the first index when a=b:otherwise 0

      integer ifind,n,i

      character a*(*),b(n)*(*)
      ifind=0
      do i=1,n
        if(a.eq.b(i)) then
            ifind=i
            return
        endif
      end do
      return
      end
