C----------------------------------------------------------------------
C
      SUBROUTINE SMOOTHXYZ (N, X, Y, Z, FWHM, NITER)
C
C  PURPOSE:
C
C     SMOOTHXYZ smooths a surface in 3-space defined by N points which
C     are treated as unstructured.  (It can also be applied to structured
C     data if the points are packed.)  The data points are overwritten.
C     One or more passes through the data may be specified.
C
C  METHOD:
C
C     One-dimensional smoothing with a Gaussian kernel is applied to each
C     of X, Y, and Z.  For each target point, the 3-space distance to each
C     other point is used as the abscissa for the Gaussian curve centered
C     at the target point with specified FWHM (full width at half maximum).
C     The X, Y, Z units should thus be comparable for distances to be
C     meaningful.  If Z(X,Y) is not a 3-space coordinate, try scaling it
C     first then unscaling upon return from the smoothing routine.
C
C     Note: FWHM is usually used for smoothing rather than the sigma that
C     is familiar from probability theory:
C
C     FWHM = sigma * sqrt (8 * ln (2)) = 2 * sigma * 1.17741...
C
C  REFERENCE:
C
C     See http://www.mrc-cbu.cam.ac.uk/Imaging/smoothing.html for an
C     outline of the Gaussian kernel method.
C
C  PERFORMANCE NOTE:
C
C     Smoothing of a surface triangulation of the Itokawa asteroid with
C     nearly 100,000 points took about 4 minutes for one iteration on a
C     2014 workstation with plenty of memory.
C
C  HISTORY:
C     09/26/00  D. Saunders  Initial implementation prompted by the
C                            above web page.
C     12/30/14     "    "    Automatic arrays for X/Y/ZSMOOTH are unwise
C                            for large N, so make them allocatable.
C
C  AUTHOR:  David Saunders, ELORET/NASA Ames Researcg Center, CA.
C                           Now with ERC, Inc. at ARC.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   N                  ! Number of points to be smoothed

      REAL, DIMENSION (N), INTENT (INOUT) ::
     >   X, Y, Z            ! Input with points to be smoothed;
     >                      ! Output with the smoothed points

      REAL, INTENT (IN) ::
     >   FWHM               ! Controls the extent of the weighted averaging
                            ! at each point; units are the distance from the
                            ! current point in (x,y,z) space; since FWHM is
                            ! about 2.35 sigma, points beyond 1.3 * FWHM have
                            ! negligible influence

      INTEGER, INTENT (IN) ::
     >   NITER              ! Number of passes to make through the data

C     Local constants:

      REAL, PARAMETER ::
     >   RLN2 = 0.69314718055994530941, ! ln (2)
     >   HALF = 0.5,
     >   ONE  = 1.,
     >   TWO  = 2.,
     >   ZERO = 0.

C     Local variables:

      INTEGER
     >   I, ITER, J
      REAL
     >   DIST, GAUSS, GSUM, RSIGMA, X0, Y0, Z0, XSUM, YSUM, ZSUM
      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   XSMOOTH, YSMOOTH, ZSMOOTH

C     Execution:

      ALLOCATE (XSMOOTH(N), YSMOOTH(N), ZSMOOTH(N))

      RSIGMA = (TWO * SQRT (TWO * RLN2)) / FWHM

      ! The factor of 1/(sigma * sqrt (2 pi)) cancels out ...

      DO ITER = 1, NITER ! For each pass

         DO I = 1, N ! Smooth each point

            GSUM = ZERO
            XSUM = ZERO
            YSUM = ZERO
            ZSUM = ZERO
            X0 = X(I)
            Y0 = Y(I)
            Z0 = Z(I)

            DO J = 1, N ! Scan all points
               DIST  = SQRT ((X(J) - X0) ** 2 + (Y(J) - Y0) ** 2 +
     >                       (Z(J) - Z0) ** 2)
               GAUSS = EXP (-HALF * (DIST * RSIGMA) ** 2)

               GSUM  = GSUM + GAUSS
               XSUM  = XSUM + GAUSS * X(J)
               YSUM  = YSUM + GAUSS * Y(J)
               ZSUM  = ZSUM + GAUSS * Z(J)
            END DO

            GSUM = ONE / GSUM
            XSMOOTH(I) = XSUM * GSUM
            YSMOOTH(I) = YSUM * GSUM
            ZSMOOTH(I) = ZSUM * GSUM

         END DO ! Smooth the next point

         DO I = 1, N
            X(I) = XSMOOTH(I)
            Y(I) = YSMOOTH(I)
            Z(I) = ZSMOOTH(I)
         END DO

      END DO ! Next iteration

      DEALLOCATE (XSMOOTH, YSMOOTH, ZSMOOTH)

      END SUBROUTINE SMOOTHXYZ
