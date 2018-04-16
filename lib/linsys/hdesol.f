C+----------------------------------------------------------------------
C
      SUBROUTINE HDESOL (MDIM, M, N1, A1, X, SSQMIN)
C
C  Purpose:
C
C     HDESOL (Householder DEcomposition and SOLution) solves one over-
C     determined system  Ax ~ b  where A is m x n with m >= n.  This is
C     the familiar linear least squares problem:  minimize the 2-norm of
C     Ax - b  with respect to x.
C
C  Notes:
C
C     1. If more than one right-hand-side b is to be used with the same
C        matrix A, use subroutines HDECOM, HSOLVE.
C     2. For the common polynomial case, see PNFIT and PNEVAL.
C     3. The direct factorization of A used here avoids the squaring of
C        the condition number that results from the alternative "normal
C        equations" approach, which solves the square system A'Ax = A'b.
C     4. No attempt is made here to handle rank-deficient cases.
C
C  Method:
C
C     A sequence of (orthogonal) Householder transformations (product =
C     Q) is used to triangularize matrix A one column at a time:
C                                                                T
C        Q   Q   ... Q Q A = R            so that           A = Q R
C         n-1 n-2     2 1
C
C     where R is a right (upper) triangular m x n matrix.
C
C     Inputting the right-hand-side b as the (n+1)th column makes for
C     convenient transformation of b during the factorization.  Then
C     the optimal solution is obtained from the upper n x n portion of
C     the transformed triangular system.
C
C  Arguments:
C
C     MDIM      is the declared row dimension of the array containing
C               matrix A.
C     M         is the number of equations in the system.  M >= N.
C     N1        is N+1 where N is the number of linear parameters X(1:N)
C               being computed.  N >= 1, so N1 >= 2.
C     A1(M,N1)  is the augmented matrix (A b), where A is the matrix
C               of coefficients in columns 1:N, and B is the RHS vector
C               in column N1 = N+1.  A is assumed to be full-rank.
C               A1 is destroyed upon return.
C     X(M)      is used for work space before being output with the
C               required linear parameters in X(1) through X(N).
C     SSQMIN    is output with the square of the  2-norm of the minimum
C               residual (minimum sum of squares).
C
C  Reference:
C
C     "Linear Least Squares Solutions by Householder Transformations,"
C     by Businger and Golub (1965).  Numer. Math. 7, pp 269-276.
C
C  Acknowledgement:
C
C     G.Golub, C.Moler, M.Saunders, Stanford University, 1972.
C
C     Jan. 1989 D.Saunders  Description clarified.
C     July 2000   "   "     Fortran 90 translation.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER, INTENT (IN)    :: MDIM, M, N1
      REAL,    INTENT (INOUT) :: A1(MDIM,N1), X(M)
      REAL,    INTENT (OUT)   :: SSQMIN

C  *  Local constants:

      REAL, PARAMETER :: ZERO = 0.

C  *  Local variables:

      INTEGER  I, J, K, N
      REAL     ALPHA, BETA, GAMMA, T1

C     Execution:

      N = N1 - 1

      DO K = 1, N1

C  *     Find the reflection which zeros A1(I,K), I = K+1, ... M.

         SSQMIN = ZERO

         IF (K < N1) THEN

            DO I = K, M
               X(I) = A1(I,K)
               SSQMIN = SSQMIN + X(I) ** 2
            END DO

         ELSE ! K = N1

            IF (M > N) THEN

               DO I = K, M
                  X(I) = A1(I,K)
                  SSQMIN = SSQMIN + X(I) ** 2
               END DO

            END IF ! SSQMIN is thus left with its required value.

            EXIT

         END IF

         ALPHA = SQRT (SSQMIN)
         T1 = X(K)
         IF (T1 < ZERO) ALPHA = -ALPHA
         T1 = T1 + ALPHA
         BETA = ALPHA * T1
         X(K) = T1

C  *     Apply the reflection to the remaining columns of A1.

         DO J = K + 1, N1

            GAMMA = ZERO
            DO I = K, M
               GAMMA = GAMMA + A1(I,J) * X(I)
            END DO
            GAMMA = GAMMA / BETA

            DO I = K, M
               A1(I,J) = A1(I,J) - GAMMA * X(I)
            END DO

         END DO

         X(K) = -ALPHA

      END DO

C  *  Solve the triangular system  Rx = <upper N part of transformed b>.

      X(N) = -A1(N,N1) / ALPHA

      IF (N > 1) THEN

         DO I = N - 1, 1, -1
            T1 = A1(I,N1)
            DO J = I + 1, N
               T1 = T1 - A1(I,J) * X(J)
            END DO
            X(I) = T1 / X(I)
         END DO

      END IF

      END SUBROUTINE HDESOL
