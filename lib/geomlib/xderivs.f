C+----------------------------------------------------------------------
C
      SUBROUTINE XDERIVS (N, XP, XPP, YP, YPP, DYDX, D2YDX2, SCALE)
C
C  PURPOSE:
C
C        XDERIVS modularizes the calculation of Y vs. X derivatives
C     from parametric derivatives for X and Y vs. T.  The formulas are:
C
C           dY/dX    =  Y' / X'
C           d2Y/dX2  =  (Y" X'  -  Y' X") / X' ** 3
C
C     The case of dX/dT ~ 0. is safeguarded with a relative test.
C     (1.E-10 * SCALE is substituted for |X'| less than this value.)
C
C        See CURV2D for the (similar) curvature calculation.
C
C  ARGUMENTS:
C
C        Obvious from the description: N >= 1; all others except SCALE
C     are REAL (N).  SCALE should be the X data range, needed to make
C     the test against division by very small numbers relative.
C
C        Each 2nd derivative is calculated before each 1st derivative,
C     so the calling program may pass the same array for both DYDX and
C     D2YDX2 if only the first derivatives are required.
C
C  HISTORY:
C
C     10/25/91   DAS   Initial implementation, for program PROFILE.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N
      REAL       XP (N), XPP (N), YP (N), YPP (N), DYDX (N), D2YDX2 (N),
     >           SCALE

C     Local constants:

      REAL       EPS
      PARAMETER (EPS = 1.E-10)

C     Local variables:

      INTEGER    I
      REAL       DXDT, EPSIL

C     Execution:

      EPSIL = EPS * SCALE
      DO 10, I = 1, N
         DXDT = XP (I)
         IF (ABS (DXDT) .LT. EPSIL) DXDT = SIGN (EPSIL, DXDT)
         D2YDX2 (I) = (YPP (I) - YP (I) * XPP (I) / DXDT) / DXDT ** 2
         DYDX (I) = YP (I) / DXDT    ! Avoid cubing DXDT to delay overflow
   10 CONTINUE

      RETURN
      END
