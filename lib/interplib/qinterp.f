C***********************************************************************
C
      SUBROUTINE QINTERP (NDATA, X, Y, NEVAL, XEVAL, YEVAL, LEFT)
C
C     QINTERP performs a form of quadratic interpolation at one or more
C     target points.  The data abscissas may be increasing or decreasing.
C     Extrapolation is permitted.  Where possible, the weighted average
C     of two interpolating quadratics is used.  The rationale is that
C     quadratics may be more theoretically correct than the cubic spline
C     that one would normally make use of, and also more trustworthy in
C     the case of extrapolation.  At least 3 data points are assumed.
C
C     07/07/00  DAS  Initial implementation, applied to aero. coefs.
C
C     Author:  David Saunders, ELORET/NASA Ames, Moffett Field, CA.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NDATA              ! Number of data points >= 3

      REAL, DIMENSION (NDATA), INTENT (IN) ::
     >   X, Y               ! Data points

      INTEGER, INTENT (IN) ::
     >   NEVAL              ! Number of evaluation points >= 1

      REAL, DIMENSION (NEVAL), INTENT (IN) ::
     >   XEVAL              ! Target abscissa(s)

      REAL, DIMENSION (NEVAL), INTENT (OUT) ::
     >   YEVAL              ! Interpolated ordinate(s)

      INTEGER, INTENT (INOUT) ::
     >   LEFT               ! Input with estimate of data interval
                            ! bracketing the first target point;
                            ! output with actual interval from the
                            ! last target point; usage is that of
                            ! the INTERVAL search utility
C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER  I, J, L
      REAL     ARROW, P, Q, XE, Y1, Y2

C     Procedure:

      EXTERNAL INTERVAL ! Search utility

C     Execution:

      L = LEFT
      ARROW = SIGN (ONE, X(2) - X(1))

      DO J = 1, NEVAL

         XE = XEVAL(J)

         CALL INTERVAL (NDATA, X, XE, ARROW, L)

         I = L

         IF (L <= 2) THEN ! Single quadratic at the low end

            I = 2

            CALL QUADRATIC (YEVAL(J), P) ! Internal procedure below

         ELSE IF (L == NDATA - 1) THEN ! Single quadratic at the high end

            CALL QUADRATIC (YEVAL(J), P)

         ELSE ! Weighted average of two quadratics

            CALL QUADRATIC (Y1, P)

            I = I + 1

            CALL QUADRATIC (Y2, Q)

            YEVAL(J) = (ONE - P) * Y1 + P * Y2

         END IF

      END DO ! Next evaluation

      LEFT = L ! May be useful if QINTERP is doing one evaluation per call


C     ******************************
C     Internal procedure for QINTERP
C     ******************************

      CONTAINS

         SUBROUTINE QUADRATIC (YE, R) ! Evaluate the quadratic through
                                      ! data points I-1, I, I+1 at X = XE
!        Argument:

         REAL, INTENT (OUT) :: YE,    ! Interpolated Y
     >                         R      ! Ratio to use in a weighted average

!        Local variables:

         REAL A, B, DEQ, DM1, DX, HEQ, HM1

!        Execution:

         DX  = XE - X(I)
         HM1 = X(I) - X(I-1)
         HEQ = X(I+1) - X(I)
         DM1 = (Y(I) - Y(I-1)) / HM1
         DEQ = (Y(I+1) - Y(I)) / HEQ
         R   = DX / HEQ

         A   = DEQ - DM1
         B   = HEQ * DM1 + HM1 * DEQ

         YE  = (A * DX + B) * (DX / (HM1 + HEQ)) + Y(I)

         END SUBROUTINE QUADRATIC

      END SUBROUTINE QINTERP
