C+---------------------------------------------------------------------
C
      SUBROUTINE FOILGRID (N, XMIN, XMAX, SWEIGHT, X)
C
C  ONE-LINER:  1-D grid distribution tailored for airfoils
C
C  DESCRIPTION:

C        FOILGRID generates N abscissas in the range [XMIN, XMAX],
C     using a composite sinusoidal + quadratic expression.  Rather
C     than including a linear component, the troublesome first few
C     increments are adjusted heuristically.  The root of the problem
C     lies in the fact that for the essentially quadratic shape function
C     near XMIN, the unadjusted dXs are of the form d, 3d, 5d, 7d, ...
C     where the leap from d to 3d is clearly too abrupt.  Something
C     more like 2d, 4d, 6d, 8d, ... = D, 2D, 3D, 4D, ... is sought
C     by scaling the initial increment and blending it in with a
C     Vinokur distribution.  Blending too many points can give significant
C     discontinuities, so the range of adjustment is varied according to
C     N, heuristically.  Checking results is recommended.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S  DIM   DESCRIPTION
C    N        I     I     -    Desired number of abscissas
C  XMIN,XMAX  R     I     -    First and last values to include
C   SWEIGHT   R     I     -    Weight applied to the sinusoidal
C                              component, in the range [0, 1].
C                              0.7 is recommended for airfoils.
C    X        R     O     N    Desired distribution of abscissas
C
C  HISTORY:
C    03/30/95   D.A.Saunders   Adaptation of James Reuther's original
C                              sinusoid + quadratic + linear hybrid, to
C                              overcome the too-small first increment
C                              while keeping inputs to a minimum.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER    N
      REAL       XMIN, XMAX, SWEIGHT, X (N)

C     Local constants:

      INTEGER    NTABLE
      REAL       HALF, ONE, PI, SCALE
      PARAMETER (NTABLE = 7,      ! Size of heuristic table for fudging
     >           HALF   = 0.5,
     >           ONE    = 1.0,
     >           PI     = 180.,   ! Make it radians if you have to
     >           SCALE  = 1.75)   ! Factor applied to initial dX

C     Local variables:

      REAL, PARAMETER :: R2D = 180.0d0 / acos(-1.0d0)
      INTEGER    I, NI
      REAL       D1, D2, RANGE, RN, SNORM, WQUAD, WSINE, XFUDGE, XNORM,
     >           TABLEN (NTABLE), TABLEX (NTABLE)
      DATA       TABLEN /   1.,  50.,  75., 100., 130., 150., 500. /,
     >           TABLEX /.030, .025, .010, .005, .003, .0025, .0024/

      SAVE       TABLEN, TABLEX  ! Established for SWEIGHT = 0.70 (recommended)

C     Execution.

C     Heuristic determination of range of X to fudge:

      RN = REAL (N)
      DO I = 1, NTABLE - 1
         IF (RN .LT. TABLEN (I + 1)) GO TO 20
      END DO

  20  XFUDGE = (TABLEX (I + 1) - TABLEX (I)) * (RN - TABLEN (I)) /
     >         (TABLEN (I + 1) - TABLEN (I)) + TABLEX (I)

      WSINE = SWEIGHT
      WQUAD = ONE - WSINE
      RANGE = XMAX - XMIN
      RN    = ONE / (RN - ONE)
      NI    = 2

      DO I = 2, N - 1
         XNORM = REAL (I - 1) * RN
         SNORM = HALF * (ONE - COS (XNORM * PI / R2D))
         XNORM = WSINE * SNORM + WQUAD * XNORM ** 2
         IF (XNORM .LT. XFUDGE) NI = I
         X (I) = XMIN + XNORM * RANGE
      END DO

C     Ensure precise end points, and adjust the first few increments
C     to overcome the too-small initial one:

      X (1) = XMIN

      IF (NI .GT. 4) THEN
         D1 = (X (2) - X (1)) * SCALE
         D2 = X (NI) - X (NI - 1)

         CALL HTDIS4 (.TRUE., X (1), X (NI), D1, D2, NI, X, -6, I)
      END IF

      X (N) = XMAX

      RETURN
      END
