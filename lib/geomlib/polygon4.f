C+----------------------------------------------------------------------
C
      SUBROUTINE POLYGON4 (X1, Y1, X2, Y2, X3, Y3, X4, Y4, XC, YC, AREA)
C
C PURPOSE: POLYGON4 calculates the centroid (XC, YC) and AREA of the
C          quadrilateral defined by four points.  AREA is positive if
C          the points are counterclockwise, else it is negative.
C
C METHOD:  The centroid coordinates are obtainable as area-weighted
C          averages from the centroids of the two triangles formed from
C          points 1, 2 and 4 and points 2, 3 and 4.  The area is the sum
C          of the areas of these two triangles.  If their sum is zero,
C          the centroid is returned as the average of the centroids of
C          the triangles.
C
C          The argument description is obvious.
C
C PROCEDURES:  POLYGON3 provides the centroid and area of each triangle.
C
C HISTORY: 09/25/92  DAS  Initial implementation.
C          09/28/93  DAS  Had to trap zero areas.
C
C AUTHOR:  D.A.Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL
     >   X1, Y1, X2, Y2, X3, Y3, X4, Y4, XC, YC, AREA
      
C     Local constants:

      REAL
     >   HALF, ZERO
      PARAMETER
     >  (HALF = 0.5E+0, ZERO = 0.E+0)

C     Local variables:

      REAL
     >   AREA1, AREA2, XC1, XC2, YC1, YC2

C     Procedures:

      EXTERNAL
     >   POLYGON3

C     Execution:

      CALL POLYGON3 (X1, Y1, X2, Y2, X4, Y4, XC1, YC1, AREA1)
      CALL POLYGON3 (X2, Y2, X3, Y3, X4, Y4, XC2, YC2, AREA2)

      AREA = AREA1 + AREA2

      IF (AREA .NE. ZERO) THEN
         XC = (XC1 * AREA1 + XC2 * AREA2) / AREA
         YC = (YC1 * AREA1 + YC2 * AREA2) / AREA
      ELSE
         XC = (XC1 + XC2) * HALF
         YC = (YC1 + YC2) * HALF
      END IF

      RETURN
      END
