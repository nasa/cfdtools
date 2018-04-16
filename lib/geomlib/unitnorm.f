C+------------------------------------------------------------------------------
C
      SUBROUTINE UNITNORM (METHOD, N, X, Y, T, XYEVAL, XNORM, YNORM)
C
C     One-liner:  Unit normals at points defining a 2-space curve
C
C     Description:
C
C        UNITNORM determines the X and Y components of the unit normals to a
C     2-space curve, at the data points defining the curve.  It makes use of
C     the gradients available from the local cubic spline routine, LCSFIT,
C     applied to X vs. T and Y vs. T.  Since results are needed at the data
C     points and nowhere in between, this is not as efficient as it might be
C     (it doesn't fully vectorize) but it avoids recoding derivatives w.r.t.
C     arc-length (possibly non-central or monotonic or for smoothly-closed
C     curves) in the presence of non-uniform data.
C
C     History:
C
C     28 Sep. 1996  D.A.Saunders  Initial implementation for application to
C                                 fuselage sections.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C ------------------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      CHARACTER METHOD * 1   ! I   'B' (loose), 'M' (monotonic) or 'C' (cyclic);
                             !     cyclic means the curve is closed smoothly and
                             !     the Nth point must match the first point
      INTEGER   N            ! I   Number of data points
      REAL      X (N), Y (N) ! I   The data points
      REAL      T (N)        ! I/S Cumulative chord lengths, normalized or not;
                             !     either T(1:N) are input or just T(1) = -1. is
                             !     input to force calculation here
      REAL      XYEVAL (N)   ! S   Work-space for spline evaluations of X & Y
      REAL      XNORM (N),   ! O   X and Y components of the unit normals; given
     >          YNORM (N)    !     applications may need to change signs of both

C-------------------------------------------------------------------------------

C     Local constants:

      REAL      ONE
      PARAMETER (ONE = 1.0E+0)

C     Local variables:

      INTEGER   I
      REAL      THETA

C     Procedures:

      EXTERNAL  CHORDS2D, LCSFIT

C     Execution:

C     Parameterize the data unless it's already been done:

      IF (T(1) .EQ. -ONE) CALL CHORDS2D (N, X, Y, .FALSE., THETA, T)

C     Set up dX/dT at the data points in XNORM.

      CALL LCSFIT (N, T, X, .TRUE., METHOD, N, T, XYEVAL, XNORM)

C     Likewise for dY/dT:

      CALL LCSFIT (N, T, Y, .TRUE., METHOD, N, T, XYEVAL, YNORM)

C     Turn the local gradients into components of unit normals:

      DO I = 1, N
         THETA = ATAN2 (YNORM (I), XNORM (I))
         XNORM (I) =  SIN (THETA)
         YNORM (I) = -COS (THETA)
      END DO

      RETURN
      END
