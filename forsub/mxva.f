      subroutine MXVA(SA,IAC,IAR,SB,IBC,SC,ICC,NRA,NCA)
      include 'mach.p'
      integer IAC,IAR,NRA,NCA,I,J,ICC,IBC
      real*8  SA(1), SB(1), SC(1)
c     INITIALIZE PRODUCT
                   DO 110 I = 1, NRA
                      SC( 1 + (I-1)*ICC ) = 0.
c           ( C(I) := 0. )
110   continue
c     MULTIPLY MATRICES FROM SA AND SB
                   DO 220 J = 1, NCA
                      DO 210 I = 1, NRA
                         SC( 1 + (I-1)*ICC )
     +                     = SC( 1 + (I-1)*ICC )
     +                       + SA( 1 + (I-1)*IAC + (J-1)*IAR )
     +                         * SB( 1 + (J-1)*IBC )
c              ( C(I) := C(I) + A(I,J)*B(J) )
210   continue
220   continue
                RETURN
                END
