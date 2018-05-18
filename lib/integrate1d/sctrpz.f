C+----------------------------------------------------------------------
C
      SUBROUTINE SCTRPZ (N, X, F, AREA)
C
C     Purpose:
C
C     SCTRPZ uses a composite trapezoidal rule to estimate the integral
C     with respect to X of F(X) between X(1) and X(N) inclusive.
C
C     Input:
C           N          Number of function values supplied. (N >= 2)
C           X(1:N)     Abscissas, assumed strictly increasing, but
C                      not necessarily equally spaced.
C           F(1:N)     Function values corresponding to X(1:N).
C
C     Output:
C           AREA       Estimate of desired integral.
C
C     Author:          David Saunders, Informatics, 10/31/79.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   X (N), F (N), AREA

C     Local constants:

      REAL
     >   HALF
      PARAMETER
     >  (HALF = 0.5E+0)

C     Local variables:

      INTEGER
     >   I

C     Execution:

      AREA = (X (2) - X (1)) * F (1) + (X (N) - X (N - 1)) * F (N)

      DO 30, I = 2, N - 1
         AREA = (X (I + 1) - X (I - 1)) * F (I) + AREA
 30   CONTINUE

      AREA = AREA * HALF

      RETURN
      END
