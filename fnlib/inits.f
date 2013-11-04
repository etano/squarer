      function inits(os,nos,eta)
      implicit real*8 (a-h,o-z)
c***begin prologue  inits
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c3a2
c***keywords  initialize,orthogonal series,special function
c***author  fullerton, w., (lanl)
c***purpose  initializes an orthogonal series so that it defines the
c            number of terms to carry in the series to meet a specified
c            error.
c***description
c
c initialize the orthogonal series so that inits is the number of terms
c needed to insure the error is no larger than eta.  ordinarily, eta
c will be chosen to be one-tenth machine precision.
c
c             input arguments --
c os     array of nos coefficients in an orthogonal series.
c nos    number of coefficients in os.
c eta    requested accuracy of series.
c***references  (none)
c***routines called  xerror
c***end prologue  inits
      dimension os(nos)
c***first executable statement  inits
c     if (nos.lt.1) call xerror ( 'inits   number of coefficients lt 1',
c    1 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
c
 20   continue
c20   if (i.eq.nos) call xerror ( 'inits   eta may be too small', 28,
c    1  1, 2)
      inits = i
c
      return
      end
