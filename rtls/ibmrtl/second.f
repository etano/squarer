      subroutine second(t)
      include 'mach.p'
      integer itime,mclock
      real*8 t
      itime=mclock()
      t=itime/100.0
      return
      end
