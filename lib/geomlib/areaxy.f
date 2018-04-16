C+----------------------------------------------------------------------
C
      SUBROUTINE AREAXY ( N, X, Y, AREA )
C
C ACRONYM:   AREA of a simple irregular polygon defined by (X,Y) pairs
C            ----                                           - -
C
C PURPOSE:   AREAXY calculates the area of the polygon defined by the
C            given coordinates.  The polygon can be any shape (even
C            reentrant) except that one part of its boundary must not
C            cross another part.
C
C METHOD:    If the polygon is not closed, then point 1 is copied as
C            point N+1.  Otherwise, points 1 and N are given identical.
C            A look at areas of triangles formed with the origin and
C            pairs of adjacent given points shows that the area is given
C            (for points 1 and N+1 equal) by
C
C               0.5 * SUM (k = 1:N) ( X(k) Y(k+1) - X(k+1) Y(k) )
C
C            Since the sign of this is positive or negative depending on
C            the direction of traversal, the absolute value is used.
C
C REFERENCE: Article in BYTE, Feb. 1987, p. 137, by Stolk and Ettershank
C            (Monash University, Clayton, Australia)
C
C ARGUMENTS:
C    ARG     DIM   TYPE I/O/S DESCRIPTION
C     N       -     I     I   Number of points supplied
C     X     N(+1)   R     I   Perimeter coordinates of polygon
C     Y       "     "     "   (with room for closing if not closed already)
C   AREA      -     R     O   Area of polygon
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C KNOWN SYSTEM DEPENDENCIES:
C   IMPLICIT NONE is not in the ANSI standard
C
C ERROR HANDLING:  None.
C
C FURTHER USAGE NOTES:
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   03/04/87   DAS    Initial implementation
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N
      REAL
     >   AREA, X (N+1), Y (N+1)

C     Local variables.

      INTEGER
     >   I, NUSED

C     Execution.

      IF ( X (1) .EQ. X (N) .AND. Y (1) .EQ. Y (N) ) THEN
         NUSED = N
      ELSE
         NUSED = N + 1
         X (NUSED) = X (1)
         Y (NUSED) = Y (1)
      END IF

      AREA = 0.E+0
      DO 20, I = 1, NUSED - 1
         AREA = AREA + X (I) * Y (I+1) - X (I+1) * Y (I)
   20 CONTINUE

      AREA = 0.5E+0 * ABS (AREA)

      RETURN
      END
