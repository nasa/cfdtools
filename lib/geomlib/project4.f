C+------------------------------------------------------------------------------
C
      SUBROUTINE PROJECT4 (X1, X2, X3, X4, XTARG, X0, DSQ, P, Q, IER)
C
C     One-liner: Orthogonal projection of a point onto a 4-sided cell in 3-space
C
C     Purpose:
C
C        PROJECT4 projects point XTARG onto the 4-sided cell defined by
C     X1, X2, X3, and X4, with a "counterclockwise" convention employed.
C     More precisely, since the quadrilateral may not be planar, the method
C     determines the (p,q) which minimizes the distance of XTARG to the
C     quadratic surface defined by
C
C        X0(p,q) = (1-p)(1-q) X1 + p(1-q) X2 + pq X3 + (1-p)q X4
C
C     where point X1 is taken as the origin for p and q.  Also returned are
C     the projected point, X0, and the (squared) distance to X0.
C
C        This version solves the nonlinear least squares problem via a
C     Newton iteration.  Starting guesses of 0.5 are used for p and q.
C
C     Usage:
C
C        Treating the quadrilateral as two triangles may be preferable (using
C     PROJECT3 twice) if the nearest cell to the target point is being sought,
C     because projection to a triangle leads to a noniterative linear method.
C     PROJECT4 allows the (p,q) from PROJECT3 to be converted to equivalent
C     fractional coordinates referred to the associated quadrilateral.  If
C     the quad. is a planar parallelogram, the two (p,q)s are the same, but
C     in general the quad. may be nonplanar, yet a best solution is provided
C     for that cell which can be adjusted for bidirectional parametric
C     interpolation within a structured surface grid.
C
C        PROJECT4 can be thought of as solving the inverse of the problem
C     solved by BILINT:  for a given target (u,v), BILINT calculates (p,q)
C     within the specified cell in 2-space (u/v space), from which (x,y,z)
C     can be interpolated, whereas here we need to determine the (p,q) or
C     (u,v) corresponding to a given (x,y,z) point near a cell in 3-space.
C     
C     History:
C
C     Nov. 1998  J. Reuther  Initial implementation of Newton method.
C      "    "    D. Saunders Gauss-Newton iteration analogous to the Newton
C                            iteration of BILINT seems OK for ~planar cells.
C     Sep. 1999    "     "   Clean-up and revised arguments after belatedly
C                            comparing the two approaches and finding that
C                            they match for planar quads., but full Newton
C                            is much more robust for nonplanar cells.
C                            Added step-halving to cover extreme cases,
C                            including reversal of uphill directions.
C     10/22/99     "     "   Spurious IER = 3 resulted if initial DP < EPS.
C                            Therefore, distinguish this from the step-halving
C                            convergence test, which also checks for DP < EPS.
C     07/25/03     "     "   Instead of reversing search direction, ensure that
C                            the Hessian is positive definite to begin with.
C                            Three-way convergence test is more thorough too.
C
C     AUTHOR: James Reuther/David Saunders, NASA Ames Research Center, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN), DIMENSION (3) ::
     >   X1, X2, X3, X4, ! (x,y,z) coordinates of the quadrilateral vertices
     >   XTARG           ! Target (x,y,z) point being projected

      REAL, INTENT (OUT) ::
     >   X0(3),          ! Coordinates of the projected point
     >   DSQ,            ! Corresponding squared distance from XTARG
     >   P, Q            ! Coefficients usable for bidirectional interpolation
                         ! as described above, referred to vertex X1

      INTEGER, INTENT (OUT) ::
     >   IER             ! 0: The iteration converged with p & q in [0,1];
                         ! 1: The iteration converged; p & q not in [0,1];
                         ! 2: The iteration limit was reached;
                         ! 3: The step halving failed - no convergence

C     Local constants:

      INTEGER, PARAMETER ::
     >   MAXITER = 40    ! Normally much less; positive definiteness needs more

      REAL, PARAMETER ::
     >   BIG = 1.E+30, HALF = 0.5, ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   N

      REAL
     >   A(2,2), DV(2), G(2), ALPHA, DELTA, DET_G, DET_MIN,
     >   DV_NORM, G_NORM, V_NORM, F0, P0, Q0, ROOT_EPS,
     >   TOL_F, TOL_F_HALF, TOL_F_3RD,
     >   AX1, AX2, AX3, AX4, DX, DX_DP, DX_DQ,
     >   AY1, AY2, AY3, AY4, DY, DY_DP, DY_DQ,
     >   AZ1, AZ2, AZ3, AZ4, DZ, DZ_DP, DZ_DQ

      LOGICAL
     >   CONVGD, REVERSED

C     Execution:

      AX1 = X1(1)
      AX2 = X2(1) - AX1
      AX3 = X4(1) - AX1
      AX4 = X3(1) - X2(1) - AX3

      AY1 = X1(2)
      AY2 = X2(2) - AY1
      AY3 = X4(2) - AY1
      AY4 = X3(2) - X2(2) - AY3

      AZ1 = X1(3)
      AZ2 = X2(3) - AZ1
      AZ3 = X4(3) - AZ1
      AZ4 = X3(3) - X2(3) - AZ3

      IF (PRECISION (TOL_F) > 10) THEN ! # decimal digits
         ROOT_EPS   = 1.E-7    ! True sqrt (eps) allows too many halvings
         TOL_F      = 1.E-12   ! Seek 12 digits in function DSQ
         TOL_F_HALF = 1.E-6
         TOL_F_3RD  = 1.E-4
      ELSE
         ROOT_EPS   = 1.E-4
         TOL_F      = 1.E-6
         TOL_F_HALF = 1.E-3
         TOL_F_3RD  = 1.E-2
      END IF

      DET_MIN = TOL_F ! Unclear what is best for det (Hessian)
      DV_NORM = BIG   ! Initialize measure of convergence and ...
      F0      = BIG   ! ... measure of reduced function value
      P       = HALF
      Q       = HALF
      IER     = 0
      ALPHA   = ZERO

C     Minimize the projected distance to the quadrilateral:

      DO N = 1, MAXITER
 
  100    CONTINUE ! Step-halving subiteration starts here

         X0(1)  = AX1 + AX2 * P + (AX3 + AX4 * P) * Q 
         X0(2)  = AY1 + AY2 * P + (AY3 + AY4 * P) * Q
         X0(3)  = AZ1 + AZ2 * P + (AZ3 + AZ4 * P) * Q

         DX     = X0(1) - XTARG(1)
         DY     = X0(2) - XTARG(2)
         DZ     = X0(3) - XTARG(3)

         DSQ    = DX ** 2 + DY ** 2 + DZ ** 2

         DX_DP  = AX2 + AX4 * Q
         DY_DP  = AY2 + AY4 * Q
         DZ_DP  = AZ2 + AZ4 * Q

         DX_DQ  = AX3 + AX4 * P
         DY_DQ  = AY3 + AY4 * P
         DZ_DQ  = AZ3 + AZ4 * P

C        Gradient vector:

         G(1)   = DX * DX_DP + DY * DY_DP + DZ * DZ_DP
         G(2)   = DX * DX_DQ + DY * DY_DQ + DZ * DZ_DQ

         G_NORM = MAX (ABS (G(1)), ABS (G(2)))
         V_NORM = MAX (ABS (P),    ABS (Q))

CCCCC    write (6, '(i3, a, 1p, e14.7, a, 0p, f9.6, a, 1p, 3e22.14)')
CCCCC>      n, '  ||g||:', g_norm, '  Step:', alpha,
CCCCC>      '  dsq, p, q:', dsq, p, q

C        Convergence tests from Practical Optimization, 8.2.3.2:

         IF (F0 - DSQ < TOL_F      * (ONE + DSQ)     .AND.
     >       G_NORM   < TOL_F_3RD  * (ONE + DSQ)     .AND.
     >       DV_NORM  < TOL_F_HALF * (ONE + V_NORM)) THEN
            EXIT
         END IF

C        Halve the step if the projected distance has not reduced:

         IF (DSQ > F0) THEN ! Not when N = 1

            IF (DV_NORM > ROOT_EPS * (ONE + V_NORM)) THEN

               ALPHA =  ALPHA * HALF
               P = P0 - ALPHA * DV(1)
               Q = Q0 - ALPHA * DV(2)
               DV_NORM = HALF * DV_NORM
               GO TO 100

            END IF

            IER = 3
            EXIT ! Step-halving failed - give up

         END IF

C        Calculate a search direction:

         A(1,1) = DX_DP ** 2 + DY_DP ** 2 + DZ_DP ** 2
         A(2,2) = DX_DQ ** 2 + DY_DQ ** 2 + DZ_DQ ** 2
         A(1,2) = DX * AX4 + DX_DP * DX_DQ +
     >            DY * AY4 + DY_DP * DY_DQ +
     >            DZ * AZ4 + DZ_DP * DZ_DQ
         A(2,1) = A(1,2) 

C        Ensure positive definiteness for a descent direction:

         DET_G  = A(1,1) * A(2,2) - A(1,2) ** 2

         IF (DET_G < DET_MIN) THEN
            DELTA = (DET_MIN - DET_G) / (A(1,1) + A(2,2))
CCCCC       write (6, '(A, 1P, 5E20.10)')
CCCCC>         '  DET_G, DELTA, A11, A22, A21: ',
CCCCC>            DET_G, DELTA, A(1,1), A(2,2), A(1,2)
            A(1,1) = A(1,1) + DELTA
            A(2,2) = A(2,2) + DELTA
         END IF

         CALL LUSOLVE (2, 2, A, G, IER) ! Solve G dv = g

         DV      = G
         DV_NORM = MAX (ABS (DV(1)), ABS (DV(2)))
         P0      = P
         Q0      = Q
         P       = P0 - DV(1)
         Q       = Q0 - DV(2)
         F0      = DSQ
         ALPHA   = ONE

      END DO


      IF (N >= MAXITER) THEN
         IER = 2            ! Too many iterations
      ELSE IF (IER /= 3) THEN
         IF (MAX (ABS (P - HALF), ABS (Q - HALF)) > HALF + ROOT_EPS)
     >      IER = 1         ! Converged to a point outside the cell
      END IF                ! Else IER = 0 from LUSOLVE or initialization

      END SUBROUTINE PROJECT4
