C+------------------------------------------------------------------------------
C
      SUBROUTINE QUINTIC2 (X1, Y1, YP1, X2, Y2, YP2, YPP2, YPPP2, COEFS)
C
C  ONE-LINER:  Calculate the quintic defined by (x,y,y') & (x,y,y',y",y"')
C
C  PURPOSE:
C
C        This variation of the earlier QUINTIC calculates the coefficients of
C     the quintic defined by two points and the 1st derivative at point 1 and
C     the first three derivatives at point 2.
C
C        The quintic is of this form, where H = X - X1:
C
C           Y(X)  =  Y1 + H (C1 + H (C2 + H (C3 + H (C4 + H C5))))
C
C        Both variations were prompted by an application for optimizing the
C     shape of the heat shield of an atmospheric entry vehicle.
C
C  HISTORY:  10/26/01  DAS  Initial implementation.
C            12/11/07   "   Variation to ensure curvature smoothness at point 2.
C
C  AUTHOR:   David Saunders, ELORET Corporation/NASA Ames Research Center, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN) ::
     >   X1, Y1, YP1, X2, Y2, YP2, YPP2, YPPP2

      REAL, INTENT (OUT) ::
     >   COEFS(5)

C     Local constants:

      INTEGER, PARAMETER ::
     >   N = 4   ! The method requires solution of a 4x4 system

C     Local variables:

      INTEGER
     >   IER

      REAL
     >   H, H2, H3, H4, A(N,N)

C     Execution:

      H  = X2 - X1
      H2 = H  * H
      H3 = H  * H2
      H4 = H2 * H2

      A(1,1) = H2
      A(2,1) = H + H
      A(3,1) = 2.
      A(4,1) = 0.
      A(1,2) = H3
      A(2,2) = H2 * 3.
      A(3,2) = H  * 6.
      A(4,2) = 6.
      A(1,3) = H4
      A(2,3) = H3 * 4.
      A(3,3) = H2 * 12.
      A(4,3) = H  * 24.
      A(1,4) = H4 * H
      A(2,4) = H4 * 5.
      A(3,4) = H3 * 20.
      A(4,4) = H2 * 6.

      COEFS(1) = YP1
      COEFS(2) = Y2   - Y1   - YP1  * H
      COEFS(3) = YP2  - YP1
      COEFS(4) = YPP2
      COEFS(5) = YPPP2

      CALL LUSOLVE (N, N, A, COEFS(2), IER)

      END SUBROUTINE QUINTIC2
