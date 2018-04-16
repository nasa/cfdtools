C+---------------------------------------------------------------------
C
      SUBROUTINE TRDIAG (A, B, C, R, S, N)
C
C     TRDIAG solves a general N*N tridiagonal system which is
C     assumed to be suitably diagonally dominant...
C
C     A(*) is subdiagonal   (1st element unused here);
C     B(*) is main diagonal;
C     C(*) is superdiagonal (nth element unused here);
C     R(*) is the right-hand-side vector;
C     S(*) is returned with the solution.
C
C     B(*) is overwritten during the solution;
C     S(*) and R(*) may be the same locations, in which case
C     the RHS  R(*) is also overwritten...
C
C     REFERENCE:  Collected algorithms from CACM: Algorithm 24.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   A (N), B (N), C (N), R (N), S (N)

C     Local variables:

      INTEGER
     >   J
      REAL
     >   W

C     Constants:

      REAL
     >   ONE
      PARAMETER
     >  (ONE = 1.E+0)

      W = ONE / B (1)
      S (1) = R (1) * W

      DO J = 1, N - 1
         B (J) = C (J) * W
         W =   ONE / (B (J + 1) - A (J + 1) * B (J))
         S (J + 1) = (R (J + 1) - A (J + 1) * S (J)) * W
      END DO

      DO J = N - 1, 1, -1
         S (J) = S (J) - B (J) * S (J + 1)
      END DO

      RETURN
      END
