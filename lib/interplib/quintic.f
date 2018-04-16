C+----------------------------------------------------------------------
C
      SUBROUTINE QUINTIC (X1, Y1, YP1, YPP1, X2, Y2, YP2, YPP2, COEFS)
C
C  ONE-LINER:  Calculate coefs. of a quintic between a (x,y,y',y") pair
C
C  PURPOSE:
C
C        QUINTIC calculates the coefficients of the quintic defined by
C     two points and corresponding 1st and 2nd derivatives.
C
C        The quintic is of this form, where H = X - X1:
C
C           Y(X)  =  Y1 + H (C1 + H (C2 + H (C3 + H (C4 + H C5))))
C
C        QUINTIC was prompted by an application involving optimization
C     of the shape of the heat shield of an atmospheric entry probe.
C
C  HISTORY:  10/26/01  DAS  Initial implementation.
C
C  AUTHOR:   David Saunders, ELORET/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN) ::
     >   X1, Y1, YP1, YPP1, X2, Y2, YP2, YPP2

      REAL, INTENT (OUT) ::
     >   COEFS(5)

C     Local constants:

      INTEGER, PARAMETER ::
     >   N = 3

C     Local variables:

      INTEGER
     >   IER

      REAL
     >   H, H2, H3, H4, A(3,3)

C     Execution:

      H  = X2 - X1
      H2 = H  * H
      H3 = H  * H2
      H4 = H2 * H2

      A(1,1) = H3
      A(2,1) = H2 * 3.
      A(3,1) = H  * 6.
      A(1,2) = H4
      A(2,2) = H3 * 4.
      A(3,2) = H2 * 12.
      A(1,3) = H4 * H
      A(2,3) = H4 * 5.
      A(3,3) = H3 * 20.

      COEFS(1) = YP1
      COEFS(2) = YPP1 * 0.5
      COEFS(3) = Y2   - Y1   - YP1  * H - COEFS(2) * H2
      COEFS(4) = YP2  - YP1  - YPP1 * H
      COEFS(5) = YPP2 - YPP1

      CALL LUSOLVE (N, N, A, COEFS(3), IER)

      END SUBROUTINE QUINTIC
