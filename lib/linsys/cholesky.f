C+----------------------------------------------------------------------
C
      SUBROUTINE CHOLESKY (N, A)
C
C PURPOSE:
C        CHOLESKY forms the Cholesky factorization A = GG' of matrix A,
C     where A is real, symmetric, and stored in lower triangular form by
C     rows.  A is assumed to be positive definite in order for the
C     factorization to be guaranteed without complex arithmetic.  Row
C     or column interchanges are also redundant in this case.  A is
C     overwritten by the factor G, which is also stored by rows.
C
C ARGUMENTS:
C   ARG  DIM     TYPE I/O/S DESCRIPTION
C    N    -       I   I     Order of the system (N >= 2)
C    A N*(N+1)/2  R   I/O   Input with lower triangle of A by rows;
C                           output with triangular factor G, by rows.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C HISTORY:
C   08/10/90  DAS  Initial implementation, with complex version in mind.
C
C AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   A (N * (N + 1) / 2)

C     Local constants:

      REAL
     >   ONE, ZERO
      PARAMETER
     >  (ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, I1, I2, J, K, L
      REAL
     >   RGII, SUM

C     Execution:

C     Elements of G are calculated by columns:

      I1 = 1
      I2 = 0
      DO 50, J = 1, N
         SUM = ZERO
         I1 = I2 + 1
         I2 = I2 + J
C***     DO 20, K = 1, J - 1       ! For A (I, J) subscripting
         DO 20, K = I1, I2 - 1
            SUM = A (K) ** 2 + SUM
   20    CONTINUE

         A (I2) = SQRT (A (I2) - SUM)
         RGII = ONE / A (I2)

         L = I2 + 1
         DO 40, I = J + 1, N
            SUM = ZERO
C***        DO 30, K = 1, J - 1
            DO 30, K = I1, I2 - 1
               SUM = A (L) * A (K) + SUM
               L = L + 1
   30       CONTINUE
            A (L) = (A (L) - SUM) * RGII
            L = L + I - J + 1           ! Not an obvious increment!
   40    CONTINUE

   50 CONTINUE

      RETURN
      END
