C+------------------------------------------------------------------------------
C
      SUBROUTINE POLY6 (X1, Y1, YP1, YPP1, X2, Y2, YP2, YPP2, YPPP2,
     >                  COEFS)
C
C  PURPOSE:
C
C        POLY6 calculates the coefficients of the 6th-degree polynomial defined
C     by two points and corresponding 1st and 2nd derivatives, along with the
C     the 3rd derivative at the second point.
C
C        The polynomial is of this form, where H = X - X1:
C
C           Y(X)  =  Y1 + H (C1 + H (C2 + H (C3 + H (C4 + H (C5 + H C6)))))
C
C        Like its QUINTIC[2] predecessors, POLY6 was prompted by an application
C     involving shape optimization of an atmospheric entry vehicle forebody.
C     For an off-center apex, it enables curvature to be continuous (if not
C     necessarily smooth) across the apex.  In order to define an appropriate
C     value for YPP1 at the apex, though, using QUINTC2 first seems unavoidable.
C     The average second derivative at the apex for pairs of defining spokes
C     through the apex should produce the desired results.
C
C  HISTORY:  09/29/08  DAS  Initial implementation, adapted from QUINTIC.
C
C  AUTHOR:   David Saunders, ELORET Corporation/NASA Ames Research Center, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN) ::
     >   X1, Y1, YP1, YPP1, X2, Y2, YP2, YPP2, YPPP2

      REAL, INTENT (OUT) ::
     >   COEFS(6)

C     Local constants:

      INTEGER, PARAMETER ::
     >   N = 4
      REAL, PARAMETER ::
     >   SIX = 6.

C     Local variables:

      INTEGER
     >   IER

      REAL
     >   H, H2, H3, H4, H5, A(N,N)

C     Execution:

      H  = X2 - X1
      H2 = H  * H
      H3 = H  * H2
      H4 = H2 * H2
      H5 = H  * H4

      A(1,1) = H3
      A(2,1) = H2 * 3.
      A(3,1) = H  * SIX
      A(4,1) = SIX
      A(1,2) = H4
      A(2,2) = H3 * 4.
      A(3,2) = H2 * 12.
      A(4,2) = H  * 24.
      A(1,3) = H4 * H
      A(2,3) = H4 * 5.
      A(3,3) = H3 * 20.
      A(4,3) = H2 * 60.
      A(1,4) = H5 * H
      A(2,4) = H5 * SIX
      A(3,4) = H4 * 30.
      A(4,4) = H3 * 120.

      COEFS(1) = YP1
      COEFS(2) = YPP1 * 0.5
      COEFS(3) = Y2   - Y1   - YP1  * H - COEFS(2) * H2
      COEFS(4) = YP2  - YP1  - YPP1 * H
      COEFS(5) = YPP2 - YPP1
      COEFS(6) = YPPP2

      CALL LUSOLVE (N, N, A, COEFS(3), IER)

      END SUBROUTINE POLY6
