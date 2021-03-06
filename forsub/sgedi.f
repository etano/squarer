      SUBROUTINE SGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      include 'mach.p'

      real*8 a,t,det,work,ten
      integer lda,n,ipvt,j,k,kp1,kb,l,i,nm1,job

C***BEGIN PROLOGUE  SGEDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1,D3A1
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=SINGLE PRECISION(SGEDI-S DGEDI-D CGEDI-C),
C             DETERMINANT,INVERSE,LINEAR ALGEBRA,MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Compute the determinant and inverse of a matrix using the
C            factors computed by SGECO or SGEFA.
C***DESCRIPTION
C
C     SGEDI computes the determinant and inverse of a matrix
C     using the factors computed by SGECO or SGEFA.
C
C     On Entry
C
C        A       REAL(LDA, N)
C                the output from SGECO or SGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from SGECO or SGEFA.
C
C        WORK    REAL(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     REAL(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if SGECO has set RCOND .GT. 0.0 or SGEFA has set
C        INFO .EQ. 0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS SAXPY,SSCAL,SSWAP
C     Fortran ABS,MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SAXPY,SSCAL,SSWAP
C***END PROLOGUE  SGEDI
      dimension A(LDA,N),DET(2),WORK(N) ,IPVT(N)
C
C     COMPUTE DETERMINANT
C
C***FIRST EXECUTABLE STATEMENT  SGEDI
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0
         DET(2) = 0.0
         TEN = 10.0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0/A(K,K)
            T = -A(K,K)
            CALL SSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0
               CALL SAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL SAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL SSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
