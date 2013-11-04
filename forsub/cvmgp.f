      function cvmgp(a,b,c)
      include 'mach.p'
      real*8 cvmgp,a,b,c
      if(c.gt.0) then
              cvmgp=a
      else
              cvmgp=b
      endif
      return
      end
