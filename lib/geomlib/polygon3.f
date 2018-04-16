C+----------------------------------------------------------------------
C
      SUBROUTINE POLYGON3 (X1, Y1, X2, Y2, X3, Y3, XC, YC, AREA)
C
C PURPOSE: POLYGON3 calculates the centroid (XC, YC) and AREA of the
C          triangle defined by three points.  AREA is positive if
C          the points are counterclockwise, else it is negative.
C          POLYGON4 uses POLYGON3 for the quadrilateral case.
C
C METHOD:  The centroid coordinates are the average of the three given
C          coordinates in each case.  The area is given by the determinant
C
C                     1 | X1  Y1  1 |
C                     - | X2  Y2  1 |
C                     2 | X3  Y3  1 |
C
C          The argument description is obvious.
C
C HISTORY: 09/25/92  DAS  Initial implementation.
C
C AUTHOR:  D.A.Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL
     >   X1, Y1, X2, Y2, X3, Y3, XC, YC, AREA
      
C     Local constants:

      REAL
     >   HALF, THIRD
      PARAMETER
     >  (HALF  = 0.5E+0,
     >   THIRD = 1.E+0 / 3.E+0)

C     Execution:

      XC = (X1 + X2 + X3) * THIRD
      YC = (Y1 + Y2 + Y3) * THIRD
      AREA = (X1 * (Y2 - Y3) + X2 * (Y3 - Y1) + X3 * (Y1 - Y2)) * HALF

      RETURN
      END
