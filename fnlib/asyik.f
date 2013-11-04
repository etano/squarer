      subroutine asyik(x,fnu,kode,flgik,ra,arg,in,y)
c***begin prologue  asyik
c***refer to  besi,besk
c
c                    asyik computes bessel functions i and k
c                  for arguments x.gt.0.0 and orders fnu.ge.35
c                  on flgik = 1 and flgik = -1 respectively.
c
c                                    input
c
c      x    - argument, x.gt.0.0e0
c      fnu  - order of first bessel function
c      kode - a parameter to indicate the scaling option
c             kode=1 returns y(i)=        i/sub(fnu+i-1)/(x), i=1,in
c                    or      y(i)=        k/sub(fnu+i-1)/(x), i=1,in
c                    on flgik = 1.0e0 or flgik = -1.0e0
c             kode=2 returns y(i)=exp(-x)*i/sub(fnu+i-1)/(x), i=1,in
c                    or      y(i)=exp( x)*k/sub(fnu+i-1)/(x), i=1,in
c                    on flgik = 1.0e0 or flgik = -1.0e0
c     flgik - selection parameter for i or k function
c             flgik =  1.0e0 gives the i function
c             flgik = -1.0e0 gives the k function
c        ra - sqrt(1.+z*z), z=x/fnu
c       arg - argument of the leading exponential
c        in - number of functions desired, in=1 or 2
c
c                                    output
c
c         y - a vector whose first in components contain the sequence
c
c                                 written by
c                                 d. e. amos
c
c     abstract
c         asyik implements the uniform asymptotic expansion of
c         the i and k bessel functions for fnu.ge.35 and real
c         x.gt.0.0e0. the forms are identical except for a change
c         in sign of some of the terms. this change in sign is
c         accomplished by means of the flag flgik = 1 or -1.
c***routines called  r1mach
c***end prologue  asyik
c
      integer in, j, jn, k, kk, kode, l
      real*8 ak,ap,arg,c, coef,con,etx,flgik,fn, fnu,gln,ra,s1,s2,
     1 t, tol, t2, x, y, z
      real*8 r1mach
      dimension y(1), c(65), con(2)
      data con(1), con(2)  /
     1        3.98942280401432678e-01,    1.25331413731550025e+00/
      data c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10),
     1     c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18),
     2     c(19), c(20), c(21), c(22), c(23), c(24)/
     3       -2.08333333333333e-01,        1.25000000000000e-01,
     4        3.34201388888889e-01,       -4.01041666666667e-01,
     5        7.03125000000000e-02,       -1.02581259645062e+00,
     6        1.84646267361111e+00,       -8.91210937500000e-01,
     7        7.32421875000000e-02,        4.66958442342625e+00,
     8       -1.12070026162230e+01,        8.78912353515625e+00,
     9       -2.36408691406250e+00,        1.12152099609375e-01,
     1       -2.82120725582002e+01,        8.46362176746007e+01,
     2       -9.18182415432400e+01,        4.25349987453885e+01,
     3       -7.36879435947963e+00,        2.27108001708984e-01,
     4        2.12570130039217e+02,       -7.65252468141182e+02,
     5        1.05999045252800e+03,       -6.99579627376133e+02/
      data c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32),
     1     c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40),
     2     c(41), c(42), c(43), c(44), c(45), c(46), c(47), c(48)/
     3        2.18190511744212e+02,       -2.64914304869516e+01,
     4        5.72501420974731e-01,       -1.91945766231841e+03,
     5        8.06172218173731e+03,       -1.35865500064341e+04,
     6        1.16553933368645e+04,       -5.30564697861340e+03,
     7        1.20090291321635e+03,       -1.08090919788395e+02,
     8        1.72772750258446e+00,        2.02042913309661e+04,
     9       -9.69805983886375e+04,        1.92547001232532e+05,
     1       -2.03400177280416e+05,        1.22200464983017e+05,
     2       -4.11926549688976e+04,        7.10951430248936e+03,
     3       -4.93915304773088e+02,        6.07404200127348e+00,
     4       -2.42919187900551e+05,        1.31176361466298e+06,
     5       -2.99801591853811e+06,        3.76327129765640e+06/
      data c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56),
     1     c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64),
     2     c(65)/
     3       -2.81356322658653e+06,        1.26836527332162e+06,
     4       -3.31645172484564e+05,        4.52187689813627e+04,
     5       -2.49983048181121e+03,        2.43805296995561e+01,
     6        3.28446985307204e+06,       -1.97068191184322e+07,
     7        5.09526024926646e+07,       -7.41051482115327e+07,
     8        6.63445122747290e+07,       -3.75671766607634e+07,
     9        1.32887671664218e+07,       -2.78561812808645e+06,
     1        3.08186404612662e+05,       -1.38860897537170e+04,
     2        1.10017140269247e+02/
c***first executable statement  asyik
      tol = r1mach(3)
      tol = max(tol,1.0e-15)
      fn = fnu
      z  = (3.0e0-flgik)/2.0e0
      kk = int(z)
      do 50 jn=1,in
        if (jn.eq.1) go to 10
        fn = fn - flgik
        z = x/fn
        ra = sqrt(1.0e0+z*z)
        gln = log((1.0e0+ra)/z)
        etx = float(kode-1)
        t = ra*(1.0e0-etx) + etx/(z+ra)
        arg = fn*(t-gln)*flgik
   10   coef = exp(arg)
        t = 1.0e0/ra
        t2 = t*t
        t = t/fn
        t = sign(t,flgik)
        s2 = 1.0e0
        ap = 1.0e0
        l = 0
        do 30 k=2,11
          l = l + 1
          s1 = c(l)
          do 20 j=2,k
            l = l + 1
            s1 = s1*t2 + c(l)
   20     continue
          ap = ap*t
          ak = ap*s1
          s2 = s2 + ak
          if (max(abs(ak),abs(ap)) .lt. tol) go to 40
   30   continue
   40   continue
      t = abs(t)
      y(jn) = s2*coef*sqrt(t)*con(kk)
   50 continue
      return
      end
