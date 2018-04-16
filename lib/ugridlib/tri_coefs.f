C+------------------------------------------------------------------------------
C
      SUBROUTINE TRI_COEFS (X1, X2, X3, XTARG, X0, DSQ, P, Q, R)
C
C     One-liner: Linear interpolation coefficients within a 3-space triangle
C
C     Purpose:
C
C        TRI_COEFS calculates coefficients for linear interpolation within a
C     triangle in 3-space.  It projects point XTARG onto the [plane of the]
C     triangle defined by points X1, X2, and X3.  It returns the foot X0 of the
C     perpendicular line from XTARG to the plane of the triangle, the squared
C     length of that line, and the interpolation coefficients (p,q,r) such that
C
C        p + q + r = 1    and   p X1 + q X2 + r X3 = X0.
C
C        Note that any of p, q, r may be outside [0, 1].  See NEAREST_TRI_POINT
C     if it is desirable to adjust the initial X0 (and p, q, r) such that no
C     extrapolation is implied.  Nonlinear interpolation within a triangulation
C     actually requires the UNadjusted coefficients returned here.
C
C     Method:
C
C        Computing the shortest distance from a point to a plane is solved as a
C     linear least squares problem via direct orthogonal factorization of the
C     associated 3 x 2 system, as opposed to solving the corresponding "normal"
C     equations.
C
C     02/04/05  DAS  Added argument r to the earlier PROJECT3, and streamlined
C                    the description, in preparation for nonlinear interpolation
C                    within a surface triangulation; improved the nomenclature.
C
C     Author:   David Saunders, ELORET/NASA Ames Research Center, Moffett Field
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN), DIMENSION (3) ::
     >   X1, X2, X3,     ! (x,y,z) coordinates of the triangle vertices
     >   XTARG           ! Target (x,y,z) point being projected

      REAL, INTENT (OUT) ::
     >   X0(3),          ! Foot of the perpendicular from XTARG to the plane
     >   DSQ,            ! Corresponding squared distance
     >   P, Q, R         ! Coefficients usable for 3-pt. linear interpolation
                         ! as indicated above; P + Q + R = 1; P = 1 at X1,
                         ! Q = 1 at X2, and R = 1 at X3; some of P, Q, R are
                         ! outside [0, 1] if X0 is not inside the triangle

C     Local variables:

      INTEGER
     >   I

      REAL
     >   A(3,3), VA(3), VB(3), X(3)

C     Execution:

C     Elimination of r from the 4 x 3 system leads to a 3 x 2:

      DO I = 1, 3
         VA(I)  = X1(I) - X3(I)
         VB(I)  = X2(I) - X3(I)
         A(I,1) = VA(I)            ! LHS of Ax ~ b (destroyed by HDESOL)
         A(I,2) = VB(I)
         A(I,3) = XTARG(I) - X3(I) ! RHS
      END DO

C     Linear least squares solution of Ax ~ b by QR factorization:

      CALL HDESOL (3, 3, 3, A, X, DSQ)

      P = X(1)
      Q = X(2)
      R = 1.0 - P - Q

      X0(:) = X3(:) + VA(:) * P + VB(:) * Q  ! Since A xbar = X0 - X3

      END SUBROUTINE TRI_COEFS
