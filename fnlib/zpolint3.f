C  This is a specialized recoding of Neville's algorithm based on the
C   POLINT routine from "Numerical Recipes", but assuming N=3, and
C   ignoring the error estimation.
C  Written by Z. Sullivan, May 2004
C  This file uses a minimal number of instructions to do 3-point fitting.
      SUBROUTINE ZPOLINT3 (XA,YA,X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION XA(3),YA(3),X,Y
      DOUBLE PRECISION C1,HO,HP,HP2,W,D1,D2,DEN

      HO=XA(1)-X
      HP=XA(2)-X
      W=YA(2)-YA(1)
      DEN=HO-HP
      DEN=W/DEN
      D1=HP*DEN
      C1=HO*DEN

      HP2=XA(3)-X
      W=YA(3)-YA(2)
      DEN=HP-HP2
      DEN=W/DEN
      D2=HP2*DEN

      W=HP*DEN-D1
      DEN=HO-HP2

      IF((X+X-XA(2)-XA(3)).GT.0D0) THEN
         Y=YA(3)+D2+HP2*W/DEN
      ELSEIF((X+X-XA(1)-XA(2)).GT.0D0) THEN
         Y=YA(2)+D1+HO*W/DEN
      ELSE
         Y=YA(1)+C1+HO*W/DEN
      ENDIF

      END

