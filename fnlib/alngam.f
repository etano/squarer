      function alngam(x)
      implicit real*8 (a-h,o-z)
c***begin prologue  alngam
c***date written   770601   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7a
c***keywords  gamma function,logarithm,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the log of the absolute value of the gamma
c            function
c***description
c
c alngam(x) computes the logarithm of the absolute value of the
c gamma function at x.
c***references  (none)
c***routines called  gamma,r1mach,r9lgmc,xerror
c***end prologue  alngam
      data sq2pil / 0.9189385332 0467274e0/
      data sqpi2l / 0.2257913526 4472743e0/
      data pi     / 3.1415926535 8979324e0/
      data xmax, dxrel / 0., 0. /
c***first executable statement  alngam
      if (xmax.ne.0.) go to 10
      xmax = r1mach(2)/log(r1mach(2))
      dxrel = sqrt (r1mach(4))
c
 10   y = abs(x)
      if (y.gt.10.0) go to 20
c
c alog (abs (gamma(x))) for  abs(x) .le. 10.0
c
      alngam = log (abs (gamma(x)))
      return
c
c alog (abs (gamma(x))) for abs(x) .gt. 10.0
c
20    continue
c20   if (y.gt.xmax) call xerror ( 'alngam  abs(x) so big alngam overflo
c    1ws', 38, 2, 2)
c
      if (x.gt.0.) alngam = sq2pil + (x-0.5)*log(x) - x + r9lgmc(y)
      if (x.gt.0.) return
c
      sinpiy = abs (sin(pi*y))
c     if (sinpiy.eq.0.) call xerror ( 'alngam  x is a negative integer',
c    1  31, 3, 2)
c
c     if (abs((x-aint(x-0.5))/x).lt.dxrel) call xerror ( 'alngam  answer
c    1 lt half precision because x too near negative integer', 68, 1, 1)
c
      alngam = sqpi2l + (x-0.5)*log(y) - x - log(sinpiy) - r9lgmc(y)
      return
c
      end
