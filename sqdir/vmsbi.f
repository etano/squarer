      subroutine vmsbi(x,bsl,eps,numerm)
      implicit none
      integer numerm,nnu,nm1,n
      real*8 bsl(numerm),rat(128),d(128),v(128),test(128)
      real*8 x,eps,rx,rx2,ot,xi,xi2,xi3,xi4,aa,emax
c bsl(i,x)=e**-x*sqrt(2pix)*i(i-.5,x) ; i is modified spherical 
c                                            bessel function
c for i ranging from 1 to numerm; eps is convergance criterion
ccc
c modifications: asymptotic breakpoint unstable for high l - adjust
c  such that critical r is outside of box.
c 4/23/96 -made double precision constants.
      if(x.le.35.*numerm) then

       if(numerm.gt.128)stop
       rx=2.d0/x
       rx2=1.d0/x
       do  nnu=2,numerm
          d(nnu)=1.d20*nnu*rx+rx2
          v(nnu)=1.d0/(nnu*rx-rx2)
          rat(nnu)=1.d0/(nnu*rx-rx2)
       enddo
c iterate until convergence
       do  nm1=1,400
        ot=nm1-.5d0
        do  nnu=2,numerm
           d(nnu)=(nnu+ot)*rx+1.d0/d(nnu)
           v(nnu)=1.d0/((nnu+ot)*rx+v(nnu))
           test(nnu)=d(nnu)*v(nnu)
           rat(nnu)=rat(nnu)*test(nnu)
        enddo

        if(nm1.gt.5) then
c test for convergence
        emax=0.d0
        do  nnu=2,numerm
           emax=max(emax,abs(abs(test(nnu))-1.d0))
        enddo
        if(emax.lt.eps) go to 15
        endif
        enddo

15     continue

       do  nnu=2,numerm
          bsl(nnu)=bsl(nnu-1)*rat(nnu)
       enddo

      else

c asymptotic expansion
       xi=1.d0/x
       xi2=xi*xi
       xi3=xi2*xi
       xi4=xi3*xi
       aa=bsl(1)
       do n=1,numerm-1
        bsl(n+1)=aa*exp(-.5d0*n*(n+1)*xi)*(1.d0-.25d0*n*(n+1)*xi2+
     $     .041666667d0*n*(n+1)*(n-2)*(n+3)*xi3+
     $     .03125d0*n*(n+1)*(5*n*(n+5)-12)*xi4)
       enddo

      endif

      return
      end
