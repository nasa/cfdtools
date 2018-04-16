!----------------------------------------------------------------------
!
      SUBROUTINE SOLVE (N, NDIM, A, B, IP)
!
!     Solves  AX = B  using the LU-decomposition from DECOMP.
!     The right-hand-side is overwritten with the solution.
!     Call SOLVE once for each RHS.
!
!     INPUT :
!           N   : Order of matrix A; N >= 2 
!           NDIM: Row dimension of A declared in the calling program
!           A   : Factorized matrix obtained from DECOMP
!           B   : Right-hand-side vector
!           IP  : Pivot vector from DECOMP;
!                 if DECOMP set IP(N) to zero, SOLVE should not be used
!
!     OUTPUT:
!           B   : Solution vector X
!
!     REFERENCE:  Forsythe, Malcolm, and Moler, 1972.
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
      REAL,    INTENT (IN)    :: A(NDIM,N)
      REAL,    INTENT (INOUT) :: B(N)
      INTEGER, INTENT (IN)    :: IP(N)

!     Local variables:

      INTEGER :: I, K, M
      REAL    :: T

!     Execution:

!     We have       LU X  =  P B.
!     First, solve  L Y   =  P B:

      DO K = 1, N - 1
         M    = IP(K)
         T    = B(M)
         B(M) = B(K)
         B(K) = T
         DO I = K + 1, N
            B(I) = B(I) + A(I,K)*T
         END DO
      END DO

!     Now solve  U X  =  Y:

      DO K = N, 2, -1
         T    = B(K) / A(K,K)
         B(K) = T
         DO I = 1, K - 1
            B(I) = B(I) - A(I,K)*T
         END DO
      END DO

      B(1) = B(1) / A(1,1)

      END SUBROUTINE SOLVE
