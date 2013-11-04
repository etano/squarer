 
      subroutine MXV(A,NRA,B,NCA,C)
      include 'mach.p'
      integer NRA,NCA,I,J
c Based on Cray scilib routine.  a must be of form a(row,col).
      real*8  A(NRA,NCA), B(NCA),C(NRA)
c     INITIALIZE PRODUCT
                   DO 110 I = 1, NRA
                      C(I) = 0.
c           ( C(I) := 0. )
110   continue
c     MULTIPLY MATRICES FROM SA AND SB
                   DO 220 J = 1, NCA
                      DO 210 I = 1, NRA
                         C(I) = C(I) +  A(I,J) * B(J)
c              ( C(I) := C(I) + A(I,J)*B(J) )
210   continue
220   continue
                RETURN
                END
