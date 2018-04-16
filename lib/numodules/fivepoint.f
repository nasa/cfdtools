C+----------------------------------------------------------------------
C
      SUBROUTINE FIVEPOINT (K, X, Y, D)
C
C     FIVEPOINT calculates 5-point finite difference estimates of the
C     1st, 2nd, 3rd, and 4th derivatives of y with respect to x at the
C     middle point, for the case of nonuniformly spaced points.
C     It was prompted by the need for 4th derivatives; an analytic
C     formula appears intractable even though the uniform case is just
C
C     y''''(k) = (y(k-2) - 4y(k-1) + 6y(k) - 4y(k+1) + y(k+2)) / h^4
C
C     Defining h(k) = x(k+1) - x(k), and g(k) = x(k+2) - x(k), we apply
C     Taylor's expansion for y(x+dx) about x(k) four times for dx =
C     h(k), -h(k-1), g(k), and -g(k-2), producing 4 linear equations in
C     the 4 unknown derivatives, which are solved numerically.
C
C     12/30/98  DAS  Initial implementation.  Is there a better way?
C
C     David Saunders, Raytheon/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: K     ! Index of central point
      REAL,    INTENT (IN)  :: X(*), ! Abscissas and ordinates; K is assumed
     >                         Y(*)  ! to be at least 2 from each end
      REAL,    INTENT (OUT) :: D(4)  ! First 4 derivatives of Y w.r.t. X at X(K)

C     Local constants:

      REAL, PARAMETER :: C1 = 24., C2 = 12., C3 = 4. ! Avoid divides

C     Local variables:

      INTEGER  IER, J, N, L(4)
      REAL     A(4,4), H

C     Storage:

      DATA     L /-2, -1, 1, 2/  ! Offsets of dX from central point

C     Execution:

      DO N = 1, 4 ! The terms in Taylor's expansion for each h
         J = K + L(N)
         H =     X(J) - X(K)
         D(N) = (Y(J) - Y(K)) * C1
         A(N,1) = H      * C1
         A(N,2) = H ** 2 * C2
         A(N,3) = H ** 3 * C3
         A(N,4) = H ** 4
      END DO

      CALL LUSOLVE (4, 4, A, D, IER)

      END SUBROUTINE FIVEPOINT
