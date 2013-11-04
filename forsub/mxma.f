 
      subroutine MXMA(SA,IAC,IAR,SB,IBC,IBR,SC,ICC,ICR,NRP,M,NCP)
      include 'mach.p'
      integer IAC,IAR,IBC,IBR,ICC,ICR,NRP,M,NCP,K,I,J
      real*8  SA(1), SB(1), SC(1)
c     INITIALIZE PRODUCT
                DO 120 K = 1, NCP
                   DO 110 I = 1, NRP
                      SC( 1 + (I-1)*ICC + (K-1)*ICR ) = 0.
c           ( C(I,K) := 0. )
110   continue
120   continue
c     MULTIPLY MATRICES FROM SA AND SB
                DO 230 K = 1, NCP
                   DO 220 J = 1, M
                      DO 210 I = 1, NRP
                         SC( 1 + (I-1)*ICC + (K-1)*ICR )
     +                     = SC( 1 + (I-1)*ICC + (K-1)*ICR )
     +                       + SA( 1 + (I-1)*IAC + (J-1)*IAR )
     +                         * SB( 1 + (J-1)*IBC + (K-1)*IBR )
c              ( C(I,K) := C(I,K) + A(I,J)*B(J,K) )
210   continue
220   continue
230   continue
                RETURN
                END
