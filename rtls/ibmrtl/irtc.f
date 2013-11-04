
      integer function irtc()
      include 'mach.p'
c     integer*4 stime, t(9)
c     call ltime(stime,t)
c     write (*,*) stime,t
c     k=0
c     do 10 l=1,9
c10    k=k+t(l)
c     irtc=k
      integer itime
      call time(itime)
      irtc=itime
      end
