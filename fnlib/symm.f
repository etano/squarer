
      subroutine symm(a,b,d)
      include 'mach.p'
      real*8 a,b,d,at
c multiplies by a by d and stores back in a and b
      at=d*a
      a=at
      b=at
      return
      end
