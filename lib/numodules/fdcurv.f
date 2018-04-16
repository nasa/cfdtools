C+----------------------------------------------------------------------
C
      SUBROUTINE FDCURV (DYDX, D2YDX2, KAPPA)
C
C  PURPOSE: FDCURV (Finite-Difference CURVature estimate) returns a
C           safe-guarded value for curvature at one point on the
C           curve Y = Y (X) given the 1st and 2nd derivatives at
C           that point, using the formula
C
C               KAPPA = D2YDX2 / (1 + DYDX ** 2) ** 3/2
C
C           The sign of KAPPA is clearly that of the 2nd derivative.
C           The derivatives could well be obtained from a spline, but
C           experience shows finite differencing can be preferable.
C
C           See modules CURV2D and CURV3D for the parametric cases.
C
C  ARGUMENTS:  Obvious from the description.  KAPPA is REAL.
C
C  HISTORY: 08/17/91  Derived FDCURV from FD12K, along with FDCNTR
C                     when it was found that FD12K did not lend
C                     itself to application to one point at a time.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      REAL       DYDX, D2YDX2, KAPPA

C     Local constants:

      REAL       DYDXMAX, ONE
      PARAMETER (DYDXMAX = 1.E+10, ONE = 1.E+0)

C     Local variables:

      REAL       TERM

C     Execution:

      TERM = ONE + (MIN (ABS (DYDX), DYDXMAX)) ** 2
      KAPPA = D2YDX2 / (TERM * SQRT (TERM))

      RETURN
      END
