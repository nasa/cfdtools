C+----------------------------------------------------------------------
C
      SUBROUTINE CHSOLVE (N, A, B)
C
C PURPOSE:
C        CHSOLVE solves the linear system  Ax = b  given  A = GG' from
C     subroutine CHOLESKY, where triangular factor G is stored by rows.
C     Right-hand-side b is overwritten by the solution x.
C
C ARGUMENTS:
C   ARG  DIM     TYPE I/O/S DESCRIPTION
C    N    -       I   I     Order of the system (N >= 2)
C    A N*(N+1)/2  R   I     Lower triangle factor G, stored by rows.
C    B    N       R   I/O   Input with RHS b;
C                           output with solution x.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C HISTORY:
C   08/13/90  DAS  Initial implementation, with complex version in mind.
C
C AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   A (N * (N + 1) / 2), B (N)

C     Local constants:

      REAL
     >   ZERO
      PARAMETER
     >  (ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, I2, J, K, L
      REAL
     >   SUM

C     Execution:

C     First, solve  Gy = b:

      L = 1
      DO 30, J = 1, N
         SUM = ZERO
         DO 20, K = 1, J - 1
            SUM = A (L) * B (K) + SUM
            L = L + 1
   20    CONTINUE
         B (J) = (B (J) - SUM) / A (L)
         L = L + 1
   30 CONTINUE

C     Now solve  G'x = y:

      I2 = L - 1          ! L = N * (N + 1) / 2
      DO 50, J = N, 1, -1
         SUM = ZERO
         L = I2
         DO 40, K = N, J + 1, -1
            SUM = A (L) * B (K) + SUM
            L = L - K + 1
   40    CONTINUE
         B (J) = (B (J) - SUM) / A (L)
         I2 = I2 - 1
   50 CONTINUE

      RETURN
      END
