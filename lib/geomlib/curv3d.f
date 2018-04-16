C+----------------------------------------------------------------------
C
      SUBROUTINE CURV3D (N, XP, XPP, YP, YPP, ZP, ZPP, KAPPA)
C
C  PURPOSE: CURV2D returns non-negative values for the curvature at the
C           points on a line in (X, Y, Z) space given the 1st and 2nd
C           derivatives with respect to arc-length at those points.
C           The relevant formula at a point is
C
C              KAPPA = | X'  *  X" | / | X'| ** 3
C
C              where  X' = (XP, YP, ZP) and X" = (XPP, YPP, ZPP),
C              and "*" represents the cross-product of these vectors.
C
C  ARGUMENTS:  Obvious from the description.  N >= 1; all others are REAL (N).
C
C  HISTORY: 08/27/91  Analogue of CURV2D for the 1-point-at-a-time case.
C           10/30/91  One point only was a bad decision.  Introduced N.
C           12/21/98  The above formula applies to ANY parameterization.
C                     Using arc-length, it becomes more simply:
C                     KAPPA = SQRT (XPP ** 2 + YPP ** 2 + ZPP ** 2)
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   XP (N), XPP (N), YP (N), YPP (N), ZP (N), ZPP (N), KAPPA (N)

C     Local variables:

      INTEGER
     >   I

C     Execution:

      DO I = 1, N
C***     KAPPA (I) = SQRT (((XP (I) * YPP (I) - XPP (I) * YP (I)) ** 2 +
C*** >                      (YP (I) * ZPP (I) - YPP (I) * ZP (I)) ** 2 +
C*** >                      (ZP (I) * XPP (I) - ZPP (I) * XP (I)) ** 2)/
C*** >                     (XP (I) ** 2 + YP (I) ** 2 + ZP (I) ** 2))

         KAPPA (I) = SQRT (XPP (I) ** 2 + YPP (I) ** 2 + ZPP (I) ** 2)
      END DO

      END SUBROUTINE CURV3D
