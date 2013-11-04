      function gamma(x)
      implicit real*8 (a-h,o-z)
c***begin prologue  gamma
c***date written   770601   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7a
c***keywords  gamma function,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the gamma function.
c***description
c
c gamma computes the gamma function at x, where x is not 0, -1, -2, ....
c gamma and x are single precision.
c***references  (none)
c***routines called  csevl,gamlim,inits,r1mach,r9lgmc,xerror
c***end prologue  gamma
      dimension gcs(23)
      data gcs   ( 1) / .0085711955 90989331e0/
      data gcs   ( 2) / .0044153813 24841007e0/
      data gcs   ( 3) / .0568504368 1599363e0/
      data gcs   ( 4) /-.0042198353 96418561e0/
      data gcs   ( 5) / .0013268081 81212460e0/
      data gcs   ( 6) /-.0001893024 529798880e0/
      data gcs   ( 7) / .0000360692 532744124e0/
      data gcs   ( 8) /-.0000060567 619044608e0/
      data gcs   ( 9) / .0000010558 295463022e0/
      data gcs   (10) /-.0000001811 967365542e0/
      data gcs   (11) / .0000000311 772496471e0/
      data gcs   (12) /-.0000000053 542196390e0/
      data gcs   (13) / .0000000009 193275519e0/
      data gcs   (14) /-.0000000001 577941280e0/
      data gcs   (15) / .0000000000 270798062e0/
      data gcs   (16) /-.0000000000 046468186e0/
      data gcs   (17) / .0000000000 007973350e0/
      data gcs   (18) /-.0000000000 001368078e0/
      data gcs   (19) / .0000000000 000234731e0/
      data gcs   (20) /-.0000000000 000040274e0/
      data gcs   (21) / .0000000000 000006910e0/
      data gcs   (22) /-.0000000000 000001185e0/
      data gcs   (23) / .0000000000 000000203e0/
      data pi /3.14159 26535 89793 24e0/
c sq2pil is log (sqrt (2.*pi) )
      data sq2pil /0.91893 85332 04672 74e0/
      data ngcs, xmin, xmax, dxrel /0, 3*0.0 /
c
c lanl dependent code removed 81.02.04
c
c***first executable statement  gamma
      if (ngcs.ne.0) go to 10
c
c ---------------------------------------------------------------------
c initialize.  find legal bounds for x, and determine the number of
c terms in the series required to attain an accuracy ten times better
c than machine precision.
c
      ngcs = inits (gcs, 23, 0.1*r1mach(3))
c
      call gamlim (xmin, xmax)
      dxrel = sqrt (r1mach(4))
c
c ---------------------------------------------------------------------
c finish initialization.  start evaluating gamma(x).
c
 10   y = abs(x)
      if (y.gt.10.0) go to 50
c
c compute gamma(x) for abs(x) .le. 10.0.  reduce interval and
c find gamma(1+y) for 0. .le. y .lt. 1. first of all.
c
      n = x
      if (x.lt.0.) n = n - 1
      y = x - float(n)
      n = n - 1
      gamma = 0.9375 + csevl(2.*y-1., gcs, ngcs)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c compute gamma(x) for x .lt. 1.
c
      n = -n
c     if (x.eq.0.) call xerror ( 'gamma   x is 0', 14, 4, 2)
c     if (x.lt.0. .and. x+float(n-2).eq.0.) call xerror (  'gamma   x is
c    1 a negative integer', 31, 4, 2)
c     if (x.lt.(-0.5) .and. abs((x-aint(x-0.5))/x).lt.dxrel) call
c    1  xerror ( 'gamma   answer lt half precision because x too near ne
c    2gative integer', 68, 1, 1)
c
      do 20 i=1,n
        gamma = gamma / (x+float(i-1))
 20   continue
      return
c
c gamma(x) for x .ge. 2.
c
 30   do 40 i=1,n
        gamma = (y+float(i))*gamma
 40   continue
      return
c
c compute gamma(x) for abs(x) .gt. 10.0.  recall y = abs(x).
c
50    continue
c50   if (x.gt.xmax) call xerror ( 'gamma   x so big gamma overflows',
c    1  32, 3, 2)
c
      gamma = 0.
c     if (x.lt.xmin) call xerror ( 'gamma   x so small gamma underflows'
c    1  , 35, 2, 1)
      if (x.lt.xmin) return
c
      gamma = exp((y-0.5)*log(y) - y + sq2pil + r9lgmc(y) )
      if (x.gt.0.) return
c
c     if (abs((x-aint(x-0.5))/x).lt.dxrel) call xerror ( 'gamma   answer
c    1 lt half precision, x too near negative integer'  , 61, 1, 1)
c
      sinpiy = sin (pi*y)
c     if (sinpiy.eq.0.) call xerror ( 'gamma   x is a negative integer',
c    1 31, 4, 2)
c
      gamma = -pi / (y*sinpiy*gamma)
c
      return
      end
