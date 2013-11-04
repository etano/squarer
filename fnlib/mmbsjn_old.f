      subroutine mmbsjn (arg,n,b,ier)
      implicit real*8 (a-h,o-z)
c   purpose             - bessel function of the first kind of
c                           nonnegative integer order for
c                           real arguments
c
c   arguments    arg    - input argument. the absolute value of arg must
c                           be less than or equal to 100000. arg must be
c                n      - input parameter specifying the number of
c                           function values to be computed.
c                b      - output vector of length n containing the
c                           computed function values. b must be typed
c                           appropriately in the calling program.
c                           b(1) will contain the computed value for
c                           order zero, b(2) will contain the computed
c                           value for order 1, b(3) for order 2, etc.
c                ier    - error parameter. (output)
c                           ier = 129 + j indicates that b(i), (i=1,j)
c                             are computed to machine precision, but
c                             precision is lost for b(i), (i=j+1,n.)
c                             see the programming notes.
c
       real*8              b(n)
c                                  first executable statement
      ier = 0
      tempa = abs(arg)
      magx = int(tempa)
!hack for large arguments
      if(tempa.gt.80.d0) then
        pi=3.1415926535d0
        prefact=sqrt(2.d0/(pi*tempa))
        do i=1,n
         b(i)=profact*cos(arg-.25d0*pi*dfloat(1+i+i))
        enddo
        return
      endif

      if(n.gt.0 .and. magx.le.100000) go to 10
c                                  error return -- arg,n is out of range
      ier = 129
   10 rsign = 1.d0
      ncalc = n
c                                  use 2-term ascending series for
c                                    small arg
      tmpa4 = tempa**4
      smallx = 1.d-14
      if(tmpa4.ge.smallx) go to 20
c                                  two-term ascending series for
c                                    small arg
      tempa = 1.d0
      tempb = -.25d0*arg*arg*rsign
      b(1) = 1.+tempb
      if(n.eq.1) go to 9005
      do 15 nn=2,n
         tempa = tempa*arg/(float(2*nn-2))
         b(nn) = tempa*(1.+tempb/(float(nn)))
   15 continue
      go to 9005
c                                  initialize the calculation of p*s
   20 nbmx = n-magx
      nn = magx+1
      plast = 1.d0
      p = (dfloat(2*nn))/tempa
c                                  calculate general significance test
      test = 2.e14
      m = 0
      if(nbmx.lt.3) go to 30
c                                  calculate p*s until nn=n-1.
c                                    check for possible overflow.
      tover = 1.d100
      nstart = magx+2
      nend = n-1
      do 25 nn=nstart,nend
         pold = plast
         plast = p
         p = (dfloat(2*nn))*plast/tempa-rsign*pold
         if(p-tover) 25, 25, 35
   25 continue
      nn = nend
c                                  calculate special significance test
c                                    for nbmx.gt.2.
c
      test = max(test,sqrt(plast*1.d14)*sqrt(2.d0*p))
c
c                                  calculate p*s until significance
c                                    test passes
   30 nn = nn+1
      pold = plast
      plast = p
      p = (dfloat(2*nn))*plast/tempa-rsign*pold
      if(p.lt.test) go to 30
      if(m.eq.1) go to 55
c                                  for j*s, a strong variant of the test
c                                    is necessary. calculate it, and
c                                    calculate p*s until this test is
c                                    passed.
      m = 1
      tempb = p/plast
      tempc = (dfloat(nn+1))/tempa
      if(tempb+1.d0/tempb.gt.2.d0*tempc)tempb=tempc+sqrt(tempc**2-1.d0)
      test = test/sqrt(tempb-1.d0/tempb)
      if(p-test) 30, 55, 55
c                                  to avoid overflow, divide p*s by
c                                    tover.  calculate p*s until
c                                    abs(p).gt.1.
   35 tover = 1.d100
      p = p/tover
      plast = plast/tover
      psave = p
      psavel = plast
      nstart = nn+1
   40 nn = nn+1
      pold = plast
      plast = p
      p = (dfloat(2*nn))*plast/tempa-rsign*pold
      if(p.le.1.d0) go to 40
      tempb = (dfloat(2*nn))/tempa
      tempc = .5d0*tempb
      tempb = plast/pold
      if(tempb+1.d0/tempb.gt.2.d0*tempc)tempb=tempc+sqrt(tempc**2-1.d0)
c
c                                  calculate backward test, and find
c                                    ncalc, the highest nn such that the
c                                    test is passed.
      test = .5d0*pold*plast*(1.d0-1.d0/tempb**2)*1.d-14
      p = plast*tover
      nn = nn-1
      nend = min0(n,nn)
      do 45 ncalc=nstart,nend
         pold = psavel
         psavel = psave
         psave = (dfloat(2*nn))*psavel/tempa-rsign*pold
         if(psave*psavel-test) 45, 45, 50
   45 continue
      ncalc = nend+1
   50 ncalc = ncalc-1
c                                  the sum b(1)+2b(3)+2b(5)... is used
c                                    to normalize. m, the coefficient of
c                                    b(nn), is initialized to 2 or 0.
   55 nn = nn+1
      m = 2*nn-4*(nn/2)
c                                  initialize the backward recursion and
c                                    the normalization sum
      tempb = 0.d0
      tempa = 1.d0/p
      sum = (float(m))*tempa
      nend = nn-n
      if(nend) 80, 70, 60
c                                  recur backward via difference
c                                    equation, calculating (but not
c                                    storing) b(nn), until nn=n.
   60 do 65 l=1,nend
         nn = nn-1
         tempc = tempb
         tempb = tempa
         tempa = ((float(2*nn))*tempb)/arg-rsign*tempc
         m = 2-m
         sum = sum+(float(m))*tempa
   65 continue
c                                  store b(nn)
   70 b(nn) = tempa
      if(n.gt.1) go to 75
c                                  n=1.  since 2*tempa is added to the
c                                    sum, tempa must be subtracted
      sum = sum-tempa
      go to 110
c                                  calculate and store b(nn-1)
   75 nn = nn-1
      b(nn) = ((float(2*nn))*tempa)/arg-rsign*tempb
      if(nn.eq.1) go to 105
      m = 2-m
      sum = sum+(dfloat(m))*b(nn)
      go to 90
c                                  nn.lt.n, so store b(nn) and set
c                                  higher orders to zero
   80 b(nn) = tempa
      nend = -nend
      do 85 l=1,nend
         itemp = nn+l
         b(itemp) = 0.0d0
   85 continue
   90 nend = nn-2
      if(nend.eq.0) go to 100
c                                  calculate via difference equation and
c                                    store b(nn), until nn=2
      do 95 l=1,nend
         nn = nn-1
         b(nn) = ((dfloat(2*nn))*b(nn+1))/arg-rsign*b(nn+2)
         m = 2-m
         sum = sum+(float(m))*b(nn)
   95 continue
c                                  calculate b(1)
  100 b(1) = 2.d0*b(2)/arg-rsign*b(3)
  105 sum = sum+b(1)
c                                  normalize--if ize=1, divide sum by
c                                    cosh(arg). divide all b(nn) by sum.
  110 continue
      do 115 nn=1,n
  115 b(nn) = b(nn)/sum
      if(ncalc.eq.n) go to 9005
      ier = 129+ncalc
 9005 return
      end
