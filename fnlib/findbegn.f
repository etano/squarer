      subroutine findbegn(iu)
      include 'mach.p'
c finds the next BEGIN record or stops
c
      integer mnprm,iu,n,ipickoff
      parameter (mnprm=12)
      character  p(mnprm)*28
1     if(ipickoff(iu,p,n,mnprm).ne.0.or.n.le.0) then
           write (*,*)' data file empty ',n
           stop
      endif
      if(p(1).eq.'BEGIN') return
      go to 1
      end
