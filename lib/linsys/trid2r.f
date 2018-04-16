C+----------------------------------------------------------------------
C
      SUBROUTINE TRID2R (N, A, B, C, R, S)
C
C     TRID2R solves a general N*N tridiagonal system, which is assumed
C     to be suitably diagonally dominant, for the TWO RIGHT-HAND-SIDES
C     case.
C
C     A(*) is the subdiagonal   (1st element unused here)
C     B(*) is the main diagonal (overwritten during the solution)
C     C(*) is the superdiagonal (Nth element unused here)
C     R(*) is the 1st right-hand-side vector (overwritten with 1st soln)
C     S(*) is the 2nd right-hand-side vector (overwritten with 2nd soln)
C
C     History: Adapted from TRDIAG for the two-RHS case.
C              David Saunders, Sterling Software, 03/16/87.
C
C-----------------------------------------------------------------------

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
      R (1) = W * R (1)
      S (1) = W * S (1)

      DO J = 2, N
         B (J - 1) = C (J - 1) * W
         W  =  ONE / (B (J) - A (J) * B (J - 1))
         R (J) = W * (R (J) - A (J) * R (J - 1))
         S (J) = W * (S (J) - A (J) * S (J - 1))
      END DO

      DO J = N - 1, 1, -1
         R (J) = R (J) - B (J) * R (J + 1)
         S (J) = S (J) - B (J) * S (J + 1)
      END DO

      RETURN
      END
