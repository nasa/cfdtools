!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE LUSOLVE (N, NDIM, A, B, IER)
!
!     ONE-LINER:  Square-system solution by LU decomposition
!
!     DESCRIPTION:
!
!        LUSOLVE solves the N*N system  A x = b  by Gaussian elimination
!     with partial pivoting.  The right-hand side  b  is overwritten by the
!     solution  x.  (Matrix  A  is also overwritten.)
!
!        For more than one right-hand side, use instead DECOMP (once) and
!     SOLVE once per RHS  b.  See also DECSLV for the single-right-hand-side
!     case:  it is slightly more efficient than LUSOLVE in its processing of
!     b  as column N+1 of the matrix, but may be less convenient.
!
!        LUSOLVE was adapted from DECSLV, itself a merger of the original
!     DECOMP and SOLVE (reference below).  Note that one version of DECOMP
!     and SOLVE provides an estimate of the matrix condition number, while
!     another version avoids this extra work for applications where cond (A)
!     is of no interest.  All of these variations vary the first index in the
!     inner loops.
!
!     INPUT:
!        N   : Order of matrix  A;  N >= 2
!        NDIM: Row dim. of array  A  declared in the calling program
!        A   : N*N array containing matrix  A
!        B   : Right-hand-side N-vector  b
!
!     OUTPUT:
!        B   : Solution vector  x
!        IER : 0 means no problem was detected;
!              1 means  A  was found to be singular
!
!     REFERENCE:  Forsythe, Malcolm & Moler, 1972.
!
!     HISTORY:    07/09/93  DAS  Revision of DECSLV to avoid having to
!                                store b in an extra column of A.
!                                IER = 0 now means no error (as it should
!                                have from the start).
!                 09/12/07   "   Slight reformatting compatible with Fortran 90.
!
!     PROGRAMMING:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN)    :: N, NDIM
      REAL,    INTENT (INOUT) :: A (NDIM, N), B (N)
      INTEGER, INTENT (OUT)   :: IER

!     Local constants:

      REAL, PARAMETER :: ONE = 1.E+0, ZERO = 0.E+0

!     Local variables:

      INTEGER :: I, J, K, M
      REAL    :: T

!     Execution:

!     Perform the LU factorization, solving L y = b as we go:

      DO K = 1, N - 1

!        Determine pivot element:

         M = K
         DO I = K + 1, N
            IF (ABS (A (I, K)) > ABS (A (M, K))) M = I
         END DO

         T = A (M, K)
         A (M, K) = A (K, K)
         A (K, K) = T
         IF (T == ZERO) GO TO 90

!        Store the multipliers:

         T = -ONE / T
         DO I = K + 1, N
            A (I, K) = A (I, K) * T
         END DO

!        Apply the multipliers to current submatrix, including RHS:

         DO J = K + 1, N
            T = A (M, J)
            A (M, J) = A (K, J)
            A (K, J) = T

            DO I = K + 1, N
               A (I, J) = A (I, J) + A (I, K) * T
            END DO

         END DO

         T = B (M)
         B (M) = B (K)
         B (K) = T

         DO I = K + 1, N
            B (I) = B (I) + A (I, K) * T
         END DO

      END DO

      IF (A (N, N) == ZERO) GO TO 90

!     Back substitution (solution of U x = y):

      DO K = N, 2, -1
         T = B (K) / A (K, K)
         B (K) = T
         DO I = 1, K - 1
            B (I) = B (I) - A (I, K) * T
         END DO
      END DO

      B (1) = B (1) / A (1, 1)

      IER = 0
      GO TO 99

 90   IER = 1  ! Matrix was singular

 99   RETURN
      END
