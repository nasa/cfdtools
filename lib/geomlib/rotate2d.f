C+----------------------------------------------------------------------
C
      SUBROUTINE ROTATE2D (N, X, Y, ANGLE, P, Q)
C
C  PURPOSE:
C     ROTATE2D rotates the given point(s) (X, Y) about the point (P, Q)
C     by the indicated angle.  Input values are overwritten.
C
C  ARGUMENTS:
C  ARG    TYPE  I/O/S   DIM   DESCRIPTION
C  N       I    I        -    Number of points to rotate
C  X,      R    I/O      N    Coordinates of point(s)
C   Y
C  ANGLE   R    I        -    Angle in degrees; positive is anticlockwise
C  P,      R    I        -    Coordinates of center of rotation
C   Q
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:  06/29/90  DAS  Initial implementation
C            08/11/94  DAS  Switched to COSD, SIND form.
C
C  AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      USE TRIGD

C     Arguments:

      INTEGER
     >   N
      REAL
     >   X (N), Y (N), ANGLE, P, Q

C     Local variables:

      INTEGER
     >   I
      REAL
     >   CI, SI, XP, YP

C     Execution:

      CI = COSD (ANGLE)
      SI = SIND (ANGLE)

      DO 100, I = 1, N
         XP = X (I) - P
         YP = Y (I) - Q
         X (I) = XP * CI - YP * SI + P
         Y (I) = XP * SI + YP * CI + Q
  100 CONTINUE

      RETURN
      END
