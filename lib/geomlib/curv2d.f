C+----------------------------------------------------------------------
C
      SUBROUTINE CURV2D (N, XP, XPP, YP, YPP, KAPPA)
C
C  PURPOSE: CURV2D returns signed values for the curvature at N points
C           on a curve in (X, Y) space given the 1st and 2nd derivatives
C           with respect to arc-length at those points.  The formulas for
C           the magnitude and sign of curvature are:
C
C                   |KAPPA| = SQRT (XPP ** 2 + YPP ** 2)
C
C              sign (KAPPA) = sign (XP * YPP - XPP * YP)
C
C           See module FDCURV for the non-parametric case.
C
C  ARGUMENTS:  Obvious from the description: N >= 1; all others are REAL (N).
C
C  HISTORY: 08/27/91  Parametric form for one point at a time, analogous
C                     to FDCURV's Y = Y (X) form, which was adapted from
C                     FD12K when the latter found not to lend itself to
C                     application to one point at a time.
C           10/25/91  One point only was a bad decision.  Introduced N.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N
      REAL       XP (N), XPP (N), YP (N), YPP (N), KAPPA (N)

C     Local variables:

      INTEGER    I

C     Execution:

      DO 10, I = 1, N
         KAPPA (I) = SIGN (SQRT (XPP (I) ** 2 + YPP (I) ** 2),
     >                     XP (I) * YPP (I) - XPP (I) * YP (I))
   10 CONTINUE

      RETURN
      END
