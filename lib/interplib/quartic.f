C+----------------------------------------------------------------------
C
      SUBROUTINE QUARTIC (X1, Y1, YP1, X2, Y2, YP2, YPP2, COEFS)
C
C  ONE-LINER:  Coefs. of a quartic defined by (x,y,y') and (x,y,y',y")
C
C  PURPOSE:
C
C        QUARTIC calculates the coefficients of the quartic defined by
C     two points, the "left-hand" 1st derivative, and the "right-hand"
C     1st and 2nd derivatives.
C
C        The quartic is of this form, where H = X - X1:
C
C           Y(X)  =  Y1 + H (C1 + H (C2 + H (C3 + H C4)))
C
C        QUARTIC was prompted by an application involving optimization
C     of the shape of the heat shield of an atmospheric entry probe.
C
C  HISTORY:  11/14/01  DAS  Adaptation of QUINTIC.
C
C  AUTHOR:   David Saunders, ELORET/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN) ::
     >   X1, Y1, YP1, X2, Y2, YP2, YPP2

      REAL, INTENT (OUT) ::
     >   COEFS(4)

C     Local constants:

      INTEGER, PARAMETER ::
     >   N = 3

C     Local variables:

      INTEGER
     >   IER

      REAL
     >   H, H2, H3, A(3,3)

C     Execution:

      H  = X2 - X1
      H2 = H  * H
      H3 = H  * H2

      A(1,1) = H2
      A(2,1) = H + H
      A(3,1) = 2.
      A(1,2) = H3
      A(2,2) = H2 * 3.
      A(3,2) = H  * 6.
      A(1,3) = H2 * H2
      A(2,3) = H3 * 4.
      A(3,3) = H2 * 12.

      COEFS(1) = YP1
      COEFS(2) = Y2   - Y1   - YP1  * H
      COEFS(3) = YP2  - YP1
      COEFS(4) = YPP2

      CALL LUSOLVE (N, N, A, COEFS(2), IER)

      END SUBROUTINE QUARTIC
