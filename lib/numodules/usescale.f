C+----------------------------------------------------------------------
C
      SUBROUTINE USESCALE (MODE, NDIM, NPTS, X, Y, Z, SCALE, SHIFT, IER)
C
C One-liner:  Normalizes or denormalizes 1, 2, or 3 dimensional data.
C
C Purpose:
C   USESCALE may be used to normalize or denormalize 1, 2 or 3D data.
C   Typical use is prior to or after application of numerical methods
C   which may be sensitive to data units.
C
C   USESCALE can be used with GETSCALE, which returns safeguarded scale
C   and shift factors for normalizing one or all of X, Y, Z to [0,1].
C
C Method:
C   Each of X[, Y[, Z]] is transformed in place as indicated below for X:
C
C   Normalize:            X <-- (X - SHIFT) / SCALE
C
C   Denormalize:          X <-- X * SCALE + SHIFT
C
C Arguments:
C    NAME    DIM    TYPE I/O/S DESCRIPTION
C   MODE      -     C*1    I   'N' to normalize or 'D' to denormalize
C   NDIM      -      I     I   Number of dimensions of the data
C   NPTS      -      I     I   Number of points to be transformed
C   X        NPTS    R    I/O  + X,Y,Z coordinates
C   Y         "      "     "   | of the input/output data
C   Z         "      "     "   + (see Notes)
C   SCALE    NDIM    R     I   Scale factor(s)
C   SHIFT     "      "     "   Shift factor(s)
C   IER       -      I     O   Error code = 0 for normal return
C                                         = -1 for bad MODE
C Notes:
C   For 1D and 2D cases, arguments Y and/or Z do not apply.  The calling
C   program should pass dummies, or pass X more than once.  The 1D case
C   is degenerate, but included for completeness.
C
C Environment:
C   VAX/VMS, FORTRAN 77 with:
C   IMPLICIT NONE
C   Names up to 8 characters
C
C Author: Michael Wong, Sterling Software, Palo Alto, CA
C
C History:
C   Oct. 88  D.A.Saunders  Initial design.
C   Nov. 88  M.D.Wong      Initial coding.
C   July 89  D.A.S.        Fixed glitches in header.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER  MODE*1
      INTEGER    NDIM, NPTS, IER
      REAL       X(NPTS), Y(NPTS), Z(NPTS), SCALE(NDIM), SHIFT(NDIM)

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      INTEGER    I
      REAL       TSCALE (3), TSHIFT (3)
      LOGICAL    NORMAL

C     Execution:

      NORMAL = MODE .EQ. 'N'

      IF (.NOT. NORMAL) THEN
         IER = -1
         IF (MODE .NE. 'D') GO TO 99
      END IF

      IER = 0

C     This short loop should save code vs. reuse of scalar temporararies...

      DO 10, I = 1, NDIM
         IF (NORMAL) THEN
            TSCALE (I) = ONE / SCALE (I)
            TSHIFT (I) = -SHIFT (I) * TSCALE (I)
         ELSE
            TSCALE (I) = SCALE (I)
            TSHIFT (I) = SHIFT (I)
         END IF
   10 CONTINUE

C     Transform the X coordinates common to all cases.

      DO 20, I = 1, NPTS
         X (I) = X (I) * TSCALE (1) + TSHIFT (1)
   20 CONTINUE

      IF (NDIM .GE. 2) THEN

C        Transform Y coordinates.

         DO 40, I = 1, NPTS
            Y (I) = Y (I) * TSCALE (2) + TSHIFT (2)
   40    CONTINUE

      END IF

      IF (NDIM .EQ. 3) THEN

C        Transform Z coordinates.

         DO 60 I = 1, NPTS
            Z (I) = Z (I) * TSCALE (3) + TSHIFT (3)
   60    CONTINUE

      END IF

   99 RETURN
      END
