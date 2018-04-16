!----------------------------------------------------------------------
!
      SUBROUTINE DECOMP (N, NDIM, A, IP)
!
!     Computes the decomposition LU = PA, where  P  is a permutation
!     matrix, by Gaussian elimination with partial pivoting.
!     The input matrix is overwritten with the factorization.
!     Use subroutine SOLVE to obtain the solution of linear systems.
!     This implementation varies the first index in all the inner loops.
!
!     LUSOLVE is equivalent and preferable for the single-RHS case.
!
!     REFERENCE:  Forsythe, Malcolm, and Moler,  1972.
!
!     INPUT:
!           N   : Order of matrix A; N >= 2
!           NDIM: Row dimension of A declared in the calling program
!           A   : Matrix to be factorized
!
!     OUTPUT:
!           A(I,J), I <= J: Upper triangular factor, U
!           A(I,J), I >  J: I - L, where L is  unit lower triangular
!           IP(K),  K <  N: Index of Kth pivot row
!           IP(N)         : (-1)**(# interchanges), or 0
!
!           IP(N) = 0  upon return menas A was found to be singular,
!                      and SOLVE should not be called in this case.
!
!     | A | can be found from  A(1,1)*A(2,2)*...*A(N,N) * IP(N).
!
!     HISTORY:
!
!        c. 1973     David Saunders, Sterling Software    (Fortran 66).
!        09/24/2007    "    "    "   Minimal conversion to Fortran 90.
!
!----------------------------------------------------------------------

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN)    :: N, NDIM
      REAL,    INTENT (INOUT) :: A(NDIM,N)
      INTEGER, INTENT (OUT)   :: IP(N)

!     Local constants:

      REAL, PARAMETER :: ONE = 1.E+0, ZERO = 0.E+0

!     Local variables:

      INTEGER :: I, J, K, M
      REAL    :: T

!     Execution:

      IP(N) = 1

      DO K = 1, N - 1

!        Determine pivot element:

         M = K
         DO I = K + 1, N
            IF (ABS (A(I,K)) > ABS (A(M,K)))  M = I
         END DO

         IP(K) = M
         T     = A(M,K)

         IF (T == ZERO) THEN
            A(N,N) = ZERO
            EXIT
         END IF

         IF (M /= K) THEN
            IP(N)  = -IP(N)
            A(M,K) = A(K,K)
            A(K,K) = T
         END IF

!        Store the multipliers:

         T = -ONE / T
         DO I = K + 1, N
            A(I,K) = A(I,K)*T
         END DO

!        Apply the multipliers to current submatrix:

         DO J = K + 1, N
            T      = A(M,J)
            A(M,J) = A(K,J)
            A(K,J) = T
            IF (T == ZERO) CYCLE
            DO I = K + 1, N
               A(I,J) = A(I,J) + A(I,K)*T
            END DO
         END DO

      END DO

      IF (A(N,N) == ZERO) THEN ! Error return: matrix is singular.
         IP(N) = 0
      END IF

      END SUBROUTINE DECOMP
