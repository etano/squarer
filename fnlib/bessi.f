      subroutine BESSI(N,X,b,eps)
      include 'mach.p'

      integer N
      real*8 X,b,eps,TOX,BIP,BI,BIM,r
      integer i,iacc,M,J
      dimension b(0:n)
      real*8 bessi0
      IF (X.EQ.0.) THEN
        b(0)=1.
        do i=1,n
           b(i)=0.
        enddo
      ELSE
      IF (N.ge.2) then
        iacc=(log10(eps))**2
        TOX=2.0/X
        BIP=0.0
        BI=1.0
        M=2*((N+INT(SQRT(FLOAT(iacc*N)))))
        DO 11 J=M,1,-1
          BIM=BIP+FLOAT(J)*TOX*BI
          BIP=BI
          BI=BIM
          IF (J.le.N) B(j)=BIP
11      CONTINUE
        endif

        b(0)=bessi0(x)
        r=b(0)/bi
        do i=1,n
        b(i)=b(i)*r
        enddo
      ENDIF
      RETURN
      END
