C+------------------------------------------------------------------------------
C
      SUBROUTINE PROJECT3 (X1, X2, X3, XTARG, X0, DSQ, P, Q)
C
C     One-liner:  Orthogonal projection of a point onto a triangle in 3-space
C
C     Purpose:
C
C        PROJECT3 projects point XTARG onto the triangle defined by points
C     X1, X2, and X3.  It returns the foot X0 of the perpendicular line from
C     XTARG to the plane of the triangle, the squared length of that line,
C     and the fractional coefficients (p,q) representing the linear combination
C     of (X1 - X2) and (X3 - X2) which equals (X0 - X2).
C
C     Method:
C
C        Computing the shortest distance from a point to a plane is solved as a
C     linear least squares problem via direct orthogonal factorization of the
C     associated 3 x 2 system, as opposed to solving the corresponding "normal"
C     equations.
C
C     Usage:
C
C        A search over many triangles can use DSQ to locate the nearest triangle
C     to a given point.  If the point is within a small tolerance of being in
C     the plane of the nearest triangle, then p and q provide 3-pt. bilinear
C     coefficients for interpolating from the triangle vertices to the point,
C     as needed for perturbing surface triangulations efficiently.  The relevant
C     3-pt. formula is:
C
C                      F(p,q) = (1 - p - q) F2 + p F1 + q F3
C
C     where all of p, q, and 1 - p - q are in [0, 1] if the point corresponding
C     to (p,q) is inside the triangle.
C
C        NOTE:  A small value of DSQ on its own does not suffice to determine
C     the nearest triangle among many.  The target point could be a long way
C     from the triangle, yet still be in or near the PLANE of the triangle.
C     All of p, q, and 1 - p - q must be in [0, 1] as well as producing a DSQ
C     below some tolerance.  For searches where the target may not be in ANY of
C     the triangles tried, the best interpolation coefficients can be determined
C     by adjusting each of p and q to be in [0,1] so that no extrapolation is
C     involved, and track the triangle that gives the shortest squared distance
C     from the target point to the ADJUSTED projection point.
C
C     Reference:
C
C        See the chapter on Orthogonal Projections in Gilbert Strang's
C     "Linear Algebra and its Applications" (Academic Press)
C
C     09/14/99  DAS  Modularization of earlier ideas first pointed out by James.
C     03/09/04   "   Clarified the description.  (Successful use as a search
C                    method can depend on the logic in the loop over triangles
C                    being searched.)
C
C     Author:   David Saunders, Mark Rimlinger, James Reuther
C               NASA Ames Research Center, Moffett Field, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN), DIMENSION (3) ::
     >   X1, X2, X3,     ! (x,y,z) coordinates of the triangle vertices
     >   XTARG           ! Target (x,y,z) point being projected

      REAL, INTENT (OUT) ::
     >   X0(3),          ! Foot of the perpendicular from XTARG
     >   DSQ,            ! Corresponding squared distance
     >   P, Q            ! Coefficients usable for 3-pt. bilinear interpolation
                         ! as described above; X2 is the origin for (p,q), while
                         ! (1,0) <-> X1 and (0,1) <-> X3

C     Local variables:

      INTEGER
     >   I

      REAL
     >   A(3,3), VA(3), VB(3), X(3)

C     Execution:

      DO I = 1, 3
         VA(I)  = X1(I) - X2(I)    ! Triangle edge vectors from X2
         VB(I)  = X3(I) - X2(I)
         A(I,1) = VA(I)            ! LHS of Ax ~ b (destroyed by HDESOL)
         A(I,2) = VB(I)
         A(I,3) = XTARG(I) - X2(I) ! RHS
      END DO

C     Linear least squares solution of Ax ~ b by QR factorization:

      CALL HDESOL (3, 3, 3, A, X, DSQ)

      P = X(1)
      Q = X(2)

      DO I = 1, 3
         X0(I) = X2(I) + P * VA(I) + Q * VB(I) ! Ax gives foot of perpendicular
      END DO

      END SUBROUTINE PROJECT3
