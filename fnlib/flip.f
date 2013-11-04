      ! flips the 2nd index of 'a' according to permutation vector 'iperm'
      ! 'b' is a temporary work area

      subroutine flip(a,b,ia,iperm,n,ifl,lstart)
      include 'mach.p'
      real*8 a,b
      integer ia,iperm,n,ifl,lstart,ls0,i
      dimension a(ia,n),b(n),iperm(n+lstart-1)

      ls0=lstart-1
      ! copy a->b
      do i=1,n
        b(i)=a(1,i)
      enddo

      if(ifl.eq.1) then
        do i=1,n
          a(1,iperm(i+ls0)-ls0)=b(i)
        enddo
      else
        do i=1,n
          a(1,i)=b(iperm(i+ls0)-ls0)
        enddo
      endif

      end

      subroutine flip4(a,b,ia,iperm,n,ifl,lstart)

      include 'mach.p'
      real*4 a,b
      integer ia,iperm,n,ifl,lstart,ls0,i
      dimension a(ia,n),b(n),iperm(n+lstart-1)

      ls0=lstart-1
      ! copy a->b
      do i=1,n
        b(i)=a(1,i)
      enddo

      if(ifl.eq.1) then
        do i=1,n
          a(1,iperm(i+ls0)-ls0)=b(i)
        enddo
      else
        do i=1,n
          a(1,i)=b(iperm(i+ls0)-ls0)
        enddo
      endif
      end
