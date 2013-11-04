      subroutine second(t)
      include 'mach.p'
      real*8 t
      real*4 etime,tarray(2)
c t will be the sum of system and user time in seconds
      t = etime(tarray)
      return
      end
