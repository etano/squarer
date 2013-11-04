      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      implicit real*8 (a-h,o-z)
      INTEGER N,NP,INDX(N)
      REAL*8 A(NP,NP),B(N)
      INTEGER I,II,J,LL
      II=0
      DO 12 I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0)THEN
            DO 11 J=II,I-1
               SUM=SUM-A(I,J)*B(J)
11          CONTINUE
         ELSE IF (SUM.NE.0.) THEN
            II=I
         ENDIF
         B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13       CONTINUE
         B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
