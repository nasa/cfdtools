C+----------------------------------------------------------------------
C
      SUBROUTINE D4COEFS (N, K, X, Y, C)
C
C     D4COEFS calculates the coefficients c(1:5) for the central finite
C     difference estimate of the 4th derivative of y with respect to x
C     at the kth point of nonuniformly spaced points (x(k), y(k)).
C     For the uniform case, the coefs. are 1, -4, 6, -4, & 1 / h^4:
C
C     y''''(k) = (y(k-2) - 4y(k-1) + 6y(k) - 4y(k+1) + y(k+2)) / h^4
C
C     but a similar formula is intractable for the nonuniform case.
C     Instead, 5 equations for the case of y = 1, x, x^2, x^3, & x^4
C     can be solved numerically.  Actually, solving a 4x4 suffices,
C     and the middle coef. is derived from c1 + c2 + c3 + c4 + c5 = 0.
C
C     For points 2 or N-1, 4-point formulas are used to permit application
C     to implicit smoothing of all interior points via solution of plain
C     pentadiagonal systems.  This requires linear least squares solution
C     of 4 equations in 3 unknowns.
C
C     12/30/98  DAS  Initial approach.
C     01/02/99   "   Provided for points 2 and N-1.
C
C     David Saunders, Raytheon/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: N,   ! Number of 1-D data points
     >                         K    ! Index of target point
      REAL,    INTENT (IN)  :: X(N),! Abscissas of data points
     >                         Y(N) ! Abscissas of data points
      REAL,    INTENT (OUT) :: C(5) ! Coefs. of y(k-2), y(k-1), y(k),
                                    ! y(k+1), & y(k+2) for y''''(k);
                                    ! C(3) applies to point K for the
                                    ! K = 2 and N-1 cases as well
C     Local constants:

      REAL, PARAMETER :: ONE = 1., R24 = 24., ZERO = 0.

C     Local variables:

      INTEGER  IER, L, M
      REAL     A(4,4), XKM

C     Execution:

      IF (K == 2) THEN ! Use y''''(2) = y''''(3)

         A(1,1:4) = ONE

         DO M = 2, 3
            L = M - 1
            A(M,1) = X(1) ** L
            A(M,2) = X(2) ** L
            A(M,3) = X(3) ** L
            A(M,4) = X(4) ** L
         END DO

         A(4,1:4) = Y(1:4)

         CALL FIVEPOINT (3, X, Y, C(2)) ! C(2:5) = y', y'', y''', y'''' at k = 3

         C(1:4) = ZERO

         CALL LUSOLVE (4, 4, A, C(2), IER)

      ELSE IF (K == N - 1) THEN

         A(1,1:4) = ONE

         DO M = 2, 3
            L = M - 1
            A(M,1) = X(K-2) ** L
            A(M,2) = X(K-1) ** L
            A(M,3) = X(K)   ** L
            A(M,4) = X(N)   ** L
         END DO

         A(4,1:4) = Y(K-2:N)

         CALL FIVEPOINT (N-2, X, Y, C) ! C(1:4) = y', y'', y''', y'''' at k = n-2

         C(1:3) = ZERO
         C(5)   = ZERO

         CALL LUSOLVE (4, 4, A, C, IER)

      ELSE ! Central 5-point case accurate for quartics

         DO M = 1, 4 
            XKM = X(K) ** M
            A(M,1) = X(K-2) ** M - XKM
            A(M,2) = X(K-1) ** M - XKM
            A(M,3) = X(K+1) ** M - XKM
            A(M,4) = X(K+2) ** M - XKM
         END DO

         C(1:3) = ZERO
         C(4)   = R24

         CALL LUSOLVE (4, 4, A, C, IER)

         C(5) = C(4)
         C(4) = C(3)
         C(3) = -(C(1) + C(2) + C(4) + C(5))

      END IF

      END SUBROUTINE D4COEFS
