C+----------------------------------------------------------------------
C
      SUBROUTINE PLFIT ( NDIM, NPTS, X, Y, Z, ARC, IER )
C
C Acronym:  Parametric Linear curve "FIT" based on "arc" length 
C
C Purpose: 
C   PLFIT computes the "arc" length along an arbitrary piecewise-linear
C   curve in 1, 2, or 3 dimensional Cartesian coordinates.  (This may also
C   be thought of as the cumulative chord length approximation to true arc
C   length along smooth curves, as commonly used for parameterizing discrete
C   datasets and estimating curvature.)
C
C   PLFIT returns an array of partial sums which can be used by PLEVAL to 
C   interpolate points on the curve. This combination is useful for curves
C   with sharp corners which cannot be followed accurately by a cubic spline.
C
C   PLFIT/PLEVAL are intended for geometric applications (where X, Y, Z have
C   similar scaling).  For the case of arbitrary 2D data, where the units
C   may differ greatly, the available CHORD utility may be the wiser choice
C   because of its careful safeguards against overflow resulting from squaring.
C
C Notes:
C   For 1D and 2D cases, arguments Y and/or Z do not apply.  The calling
C   program should pass dummies, or pass X more than once.  The 1D case
C   is degenerate, but included for completeness.
C 
C Method:
C   The usual approximation to arc length is returned for each point,
C   namely cumulative chord length from the first point.  There are no
C   safeguards against overflow from squaring distances (see Purpose).
C
C Arguments:
C    NAME    DIM    TYPE I/O/S DESCRIPTION
C   NDIM      -      I     I   Number of dimensions of the curve
C   NPT       -      I     I   Number of points on the curve
C   X        NX      R     I   + X,Y,Z coordinates
C   Y         "      "     "   | of the input curve
C   Z         "      "     "   + (see Notes above)
C   ARC      NX      R     O   Partial arc lengths at each point 
C   IER       -      I     O   Error code  =0 for successful completion
C                                          -1 bad number of dimensions
C
C System Dependencies:
C   IMPLICIT NONE is non-standard.
C
C Environment:  VAX/VMS, FORTRAN 77
C
C Author: Ron Langhi, Sterling Software, Palo Alto, CA
C 
C Modification history:
C   Dec 86  R.Langhi    Original design and coding
C   Jan 87  D.Serafini  Minor mods to allow for more general usage
C   Sep 88  D.Saunders  Clarified the purpose description
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments.

      INTEGER NDIM, NPTS, IER
      REAL    X(NPTS), Y(NPTS), Z(NPTS), ARC(NPTS)

C ... Local variables.

      INTEGER I

C ... Execution.

      IER = 0
      ARC(1) = 0.E+0

C ... ARC(I) is the linear distance along the curve from the beginning
C     of the curve to point(I).

      IF ( NDIM .EQ. 3 ) THEN
         DO 20, I = 2, NPTS
            ARC(I) = ARC(I-1)
     &             + SQRT(   ( X(I) - X(I-1) )**2
     &                     + ( Y(I) - Y(I-1) )**2
     &                     + ( Z(I) - Z(I-1) )**2
     &                   )
  20     CONTINUE

      ELSE IF ( NDIM .EQ. 2 ) THEN
         DO 40, I = 2, NPTS
            ARC(I) = ARC(I-1)
     &             + SQRT(   ( X(I) - X(I-1) )**2
     &                     + ( Y(I) - Y(I-1) )**2
     &                   )
  40     CONTINUE

      ELSE IF ( NDIM .EQ. 1 ) THEN
         DO 60, I = 2, NPTS
            ARC(I) = ARC(I-1) + ABS( X(I) - X(I-1) )
  60     CONTINUE

      ELSE
         IER = -1

      END IF

      RETURN
      END
