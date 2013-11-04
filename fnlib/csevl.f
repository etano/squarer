      function csevl(x,cs,n)
      implicit real*8 (a-h,o-z)
c***begin prologue  csevl
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c3a2
c***keywords  chebyshev,fnlib,special function
c***author  fullerton, w., (lanl)
c***purpose  evaluate the n-term chebyshev series cs at x.
c***description
c
c evaluate the n-term chebyshev series cs at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973). also see fox
c and parker, chebyshev polynomials in numerical analysis, oxford press,
c page 56.
c
c       input arguments --
c x    value at which the series is to be evaluated.
c cs   array of n terms of a chebyshev series.  in eval-
c      uating cs, only half the first coefficient is summed.
c n    number of terms in array cs.
c***references  (none)
c***routines called  xerror
c***end prologue  csevl
c
       dimension cs(1)
c***first executable statement  csevl
c      if(n.lt.1) call xerror( 'csevl   number of terms le 0', 28, 2,2)
c      if(n.gt.1000) call xerror ( 'csevl   number of terms gt 1000',
c    1   31,3,2)
c      if (x.lt. -1.0 .or. x.gt. 1.0) call xerror( 'csevl   x outside (-
c    11,+1)', 25, 1, 1)
c
       b1=0.
       b0=0.
       twox=2.*x
       do 10 i=1,n
       b2=b1
       b1=b0
       ni=n+1-i
       b0=twox*b1-b2+cs(ni)
 10    continue
c
       csevl = 0.5 * (b0-b2)
c
       return
      end
