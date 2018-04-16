C+----------------------------------------------------------------------
C
      SUBROUTINE NEARPT ( X0, Y0, NXY, X, Y, NSKIP, ISKIP, INEAR )
C
C PURPOSE:  NEARPT  looks among the given NXY points (X,Y) for the point
C           nearest to given point (X0,Y0), ignoring some or none of the
C           points according to the list of indices in ISKIP(*).    This
C           capability was prompted by a requirement to  identify  three
C           noncollinear points nearest to a given point, for 2-D inter-
C           polation purposes.
C
C METHOD:   Brute force - no way of avoiding IFs inside nested loops.
C
C ARGUMENTS:
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C   X0, Y0      -       R     I   Target point
C   NXY         -       I     I   No. of points to check
C   X, Y       NXY      R     I   Points in a plane to check
C   NSKIP       -       I     I   No. of points to avoid. NSKIP >= 0.
C   ISKIP    <NSKIP>    I     I   Indices of X(*), Y(*) to avoid
C   INEAR       -       I     O   Index of X(*), Y(*) point nearest to
C                                 (X0,Y0).  INEAR=0 if NSKIP=NXY.   
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   02/13/86   DAS    Initial design and code
C
C AUTHOR: David Saunders, Sterling Software/Informatics, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Argument declarations:

      INTEGER  NXY, NSKIP, ISKIP(*), INEAR
      REAL     X0, Y0, X(NXY), Y(NXY)

C ... Local variable/array declarations:

      INTEGER  I, J
      REAL     D, DMIN

C ... Executable statements:

      IF ( NSKIP .GE. NXY ) THEN
         INEAR = 0
         RETURN
      END IF

      DMIN = 1.E+36

      DO 300 I = 1, NXY

         DO 200 J = 1, NSKIP
            IF ( I .EQ. ISKIP(J) ) GO TO 300
  200    CONTINUE

         D = ( X(I) - X0 ) ** 2  +  ( Y(I) - Y0 ) ** 2

         IF ( D .LT. DMIN ) THEN
            DMIN = D
            INEAR = I
         END IF

  300 CONTINUE

      RETURN
      END
