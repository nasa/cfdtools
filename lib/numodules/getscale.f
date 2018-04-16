C+----------------------------------------------------------------------
C
      SUBROUTINE GETSCALE (MODE, NDIM, NPTS, X, Y, Z, SCALE, SHIFT, IER)
C
C One-liner:  Determines data normalization parameters.
C
C Purpose:
C   GETSCALE calculates the scale and shift factors pertinent to normalizing
C   1, 2, or 3 dimensional data to a region near the origin ([0,1]) as may
C   be needed prior to application of certain numerical methods which are
C   sensitive to data units.  It allows the user two options:  to normalize
C   geometric data in a way that preserves the original shape, or to normalize
C   coordinates in an independent manner.  It also handles degenerate cases.
C
C   This routine can be used in conjunction with USESCALE to normalize or
C   denormalize data.
C
C Method:
C   GETSCALE finds minimum and maximum data values from input to determine
C   scale and shift factors.  Shift factors transform the minimum data value(s)
C   to the origin.  Mode 'G' produces a single scale factor for all coordinates,
C   namely the largest coordinate range.  Mode 'N' (or 'I') produces individual
C   coordinate scale factors.  Degenerate cases return unit scale factors where
C   appropriate.
C
C   In all cases, the returned parameters should be used as follows to
C   normalize:
C                 X <-- (X - SHIFT) / SCALE
C
C   and as follows to denormalize:
C
C                 X <-- X * SCALE + SHIFT
C
C Arguments:
C    NAME    DIM    TYPE I/O/S DESCRIPTION
C   MODE      -      C*1   I   'G' gives "geometric" normalization (preserving
C                              the shape).
C                              'I' or 'N' gives "independent" or "normal"
C                              normalization.
C   NDIM      -      I     I   Number of dimensions of the data.  NDIM = 1,
C                              2 or 3.
C   NPTS      -      I     I   Number of data points.  NPTS >= 1.
C   X        NPTS    R     I   + X,Y,Z coordinates
C   Y         "      "     "   | of the input data
C   Z         "      "     "   + (see Notes)
C   SCALE    NDIM    R     O   Scaling factor(s)
C   SHIFT     "      "     "   Shifting factor(s)
C   IER       -      I     O   Error code  =0 for successful completion
C                                          -1 bad number of dimensions or
C                                             invalid MODE
C
C Notes:
C   For 1D and 2D cases, arguments Y and/or Z do not apply.  The calling
C   program should pass dummies, or pass X more than once.  The 1D case is
C   degenerate, but included for completeness.
C
C External References:
C   BOUNDS   Determines data ranges
C
C Environment:
C   VAX/VMS, FORTRAN 77, with:
C   IMPLICIT NONE
C   Names up to 8 characters
C
C Author: Michael Wong, Sterling Software, Palo Alto, CA
C
C History:
C   Oct. 88  D.A.Saunders  Initial design.
C   Nov. 88  M.D.Wong      Initial coding.
C   July 88  D.A.S.        Cosmetic refinements.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER MODE*1
      INTEGER   NDIM, NPTS, IER
      REAL      X(NPTS), Y(NPTS), Z(NPTS), SCALE(NDIM), SHIFT(NDIM)

C     Local constants:

      REAL      ONE, ZERO
      PARAMETER (ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables:

      REAL      XDIFF, YDIFF, ZDIFF,
     >          XMAX, YMAX, ZMAX,
     >          XMIN, YMIN, ZMIN,
     >          XSCALE, YSCALE, ZSCALE,
     >          XSHIFT, YSHIFT, ZSHIFT, MXDIFF
      LOGICAL   GEOMETRIC

C     Execution:

      GEOMETRIC = MODE .EQ. 'G'

      IF (.NOT. GEOMETRIC) THEN
         IER = -1
         IF (MODE .NE. 'I' .AND. MODE .NE. 'N')  GO TO 999
      END IF

      IER = 0

C     Start by calculating the X parameters common to all dimensions.

C     Find minimum and maximum X values.

      XMIN = X (1)
      XMAX = XMIN

      CALL BOUNDS (NPTS, 1, NPTS, X, XMIN, XMAX)

C     Find X shift factor.

      XSHIFT = XMIN

C     Find X coordinate range.

      XDIFF = XMAX - XMIN

      IF (NDIM .EQ. 3) THEN

C        Calculate parameters for Y and Z coordinates.
C        First, find minimum and maximum Y and Z values.

         YMIN = Y(1)
         YMAX = YMIN
         ZMIN = Z(1)
         ZMAX = ZMIN

         CALL BOUNDS (NPTS, 1, NPTS, Y, YMIN, YMAX)
         CALL BOUNDS (NPTS, 1, NPTS, Z, ZMIN, ZMAX)

C        Find Y and Z shift factors.

         YSHIFT = YMIN
         ZSHIFT = ZMIN

C        Find Y and Z coordinate ranges.

         YDIFF = YMAX - YMIN
         ZDIFF = ZMAX - ZMIN
         MXDIFF = MAX (XDIFF, YDIFF, ZDIFF)

         IF (MXDIFF .EQ. ZERO) THEN

C           Handle 3D degenerate case.

            XSCALE = ONE
            YSCALE = ONE
            ZSCALE = ONE

         ELSE

C           Find X, Y and Z scaling factors.

            IF (GEOMETRIC) THEN

C              Apply geometric coordinate scaling to preserve original shape.

               XSCALE = MXDIFF
               YSCALE = XSCALE
               ZSCALE = XSCALE

            ELSE

C              Apply independent coordinate scaling.

               IF (XDIFF .EQ. ZERO) THEN
                  XSCALE = ONE
               ELSE
                  XSCALE = XDIFF
               END IF

               IF (YDIFF .EQ. ZERO) THEN
                  YSCALE = ONE
               ELSE
                  YSCALE = YDIFF
               END IF

               IF (ZDIFF .EQ. ZERO) THEN
                  ZSCALE = ONE
               ELSE
                  ZSCALE = ZDIFF
               END IF

            END IF

         END IF

C        Pack X, Y and Z shift and scale factors into output arguments.

         SHIFT (1) = XSHIFT
         SHIFT (2) = YSHIFT
         SHIFT (3) = ZSHIFT
         SCALE (1) = XSCALE
         SCALE (2) = YSCALE
         SCALE (3) = ZSCALE

      ELSE IF (NDIM .EQ. 2) THEN

C        Determine Y coordinate boundary values.

         YMIN = Y (1)
         YMAX = YMIN

         CALL BOUNDS (NPTS, 1, NPTS, Y, YMIN, YMAX)

         YSHIFT = YMIN
         YDIFF = YMAX - YMIN
         MXDIFF = MAX (XDIFF, YDIFF)

C        Handle 2D degenerate case.

         IF (MXDIFF .EQ. ZERO) THEN

            XSCALE = ONE
            YSCALE = ONE

         ELSE IF (GEOMETRIC) THEN

            XSCALE = MXDIFF
            YSCALE = XSCALE

         ELSE

            IF (XDIFF .EQ. ZERO) THEN
               XSCALE = ONE
            ELSE
               XSCALE = XDIFF
            END IF

            IF (YDIFF .EQ. ZERO) THEN
               YSCALE = ONE
            ELSE
               YSCALE = YDIFF
            END IF

         END IF

C        Pack X and Y shift and scale factors into output arguments.

         SHIFT (1) = XSHIFT
         SHIFT (2) = YSHIFT
         SCALE (1) = XSCALE
         SCALE (2) = YSCALE

      ELSE IF (NDIM .EQ. 1) THEN

C        Find X scale factor.

         IF (XDIFF .EQ. ZERO) THEN

C           Handle 1D degenerate case.

            XSCALE = ONE
         ELSE
            XSCALE = XDIFF
         END IF

C        Pack X shift and scale factors into output arguments.

         SHIFT (1) = XSHIFT
         SCALE (1) = XSCALE

      ELSE

         IER = -1  ! Invalid NDIM

      END IF

999   RETURN
      END
