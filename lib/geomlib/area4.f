C+-----------------------------------------------------------------------
C
      REAL FUNCTION AREA4 (X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3,
     >                     X4, Y4, Z4)
C
C ONE-LINER:  Area of a quadrilateral in 3-space.
C
C PURPOSE:
C     AREA4 calculates the area of a quadrilateral defined by four coplanar
C     points P1, P2, P3, and P4 in space.  The points are assumed to be in
C     either clockwise or counterclockwise order around the perimeter.
C
C ARGUMENTS:
C     ARG       TYPE    I/O      DESCRIPTION
C     X1        R      I       X coordinate of P1
C     Y1        R      I       Y     "         "
C     Z1        R      I       Z     "         "
C     .                .       .
C     .                .       .
C     .                .       .
C     Z4        R      I       Z coordinate of P4
C
C METHOD:
C     The function calculates the area of the quadrilateral as the sum of the
C     areas of the two triangles defined by diagonal P1 and P3.  The areas of
C     the triangles are calculated using the cross product method. (Reference:
C     "Standard Mathematical Tables", Beyer, P. 204.)
C
C HISTORY:
C     3/31/89     M.Wong     Initial design and implementation.
C
C AUTHOR:   Michael D. Wong, Sterling Software, Palo Alto, CA
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL       X1, X2, X3, X4,
     >           Y1, Y2, Y3, Y4,
     >           Z1, Z2, Z3, Z4

C     Local variables:

      REAL       AREASUB1, AREASUB2

C     Execution:

      AREASUB1 = SQRT ((Y1 * (Z2 - Z3) - Z1 * (Y2 - Y3) +
     >                 (Y2 * Z3 - Y3 * Z2)) ** 2
     >                               +
     >                 (Z1 * (X2 - X3) - X1 * (Z2 - Z3) +
     >                 (Z2 * X3 - Z3 * X2)) ** 2
     >                               +
     >                 (X1 * (Y2 - Y3) - Y1 * (X2 - X3) +
     >                 (X2 * Y3 - X3 * Y2)) ** 2)

      AREASUB2 = SQRT ((Y1 * (Z4 - Z3) - Z1 * (Y4 - Y3) +
     >                 (Y4 * Z3 - Y3 * Z4)) ** 2
     >                               +
     >                 (Z1 * (X4 - X3) - X1 * (Z4 - Z3) +
     >                 (Z4 * X3 - Z3 * X4)) ** 2
     >                               +
     >                 (X1 * (Y4 - Y3) - Y1 * (X4 - X3) +
     >                 (X4 * Y3 - X3 * Y4)) ** 2)

      AREA4 = 0.5 * (AREASUB1 + AREASUB2)

      END FUNCTION AREA4
