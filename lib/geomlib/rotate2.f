C+----------------------------------------------------------------------
C
      SUBROUTINE ROTATE2 (N, X, Y, XN, YN, ANGLE, P, Q)
C
C  PURPOSE:
C     ROTATE2 rotates the given point(s) (X, Y) about the point (P, Q)
C     by the indicated angle.  Input values may be overwritten.
C
C  ARGUMENTS:
C  ARG    TYPE  I/O/S   DIM   DESCRIPTION
C  N       I    I        -    Number of points to rotate
C  X,      R    I        N    Input coordinates of point(s)
C   Y
C  XN,     R      O      N    Output coordinates
C   YN
C  ANGLE   R    I        -    Angle in degrees; positive is anticlockwise
C  P,      R    I        -    Coordinates of center of rotation
C   Q
C
C  HISTORY: 06/29/90    DAS   Initial implementation (ROTATE2D).
C              ?    J.Reuther ROTATE2: Expect radians, not degrees;
C                             output locations may differ from input.
C           07/28/94    DAS   Cleaned up and renamed as ROTATE2.
C           08/15/94    DAS   Reverted to degrees: COSD, SIND, ATAN2D
C                             are now virtually standard, and preferred.
C
C  AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

C     Arguments:

      INTEGER
     >   N
      REAL
     >   X (N), Y (N), XN (N), YN (N), ANGLE, P, Q

C     Local variables:

      INTEGER
     >   I
      REAL
     >   CI, SI, XP, YP

C     Execution:

      CI = COSD (ANGLE)
      SI = SIND (ANGLE)

      DO 100, I = 1, N
         XP     = X (I) - P
         YP     = Y (I) - Q
         XN (I) = XP * CI - YP * SI + P
         YN (I) = XP * SI + YP * CI + Q
  100 CONTINUE

      RETURN
      END
