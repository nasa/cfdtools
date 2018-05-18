C+------------------------------------------------------------------------------
C
      SUBROUTINE LCSAREAS (NDATA, X, Y, METHOD, OFFSET, AREAS)
C
C  One-liner:  Local cubic spline integration of discrete data (partial sums)
C
C  Description:
C
C        LCSAREAS calculates the definite integrals up to each knot for the
C     local cubic spline defined by discrete data points in 2-space, with an
C     optional constant added to all partial integrals.  Local cubic splines
C     use 4-pt. methods to determine spline coefficients on the fly as needed.
C
C        This variation of LCSQUAD returns all partial integrals on the
C     intervals [X(1), X(I)], I = 1 : NDATA, whereas LCSQUAD returns a single
C     definite integral on the (arbitrary) interval [XA, XB].  The "CYCLIC"
C     option is also omitted here for added efficiency.
C
C  Arguments:
C
C     Name    Dimension  I/O/S  Description
C
C     NDATA              I      Length of X, Y input data arrays.
C
C     X,      (NDATA)    I      Input data coordinates.  The Xs
C     Y                         must be distinct and monotonic,
C                               either ascending or descending.
C                               (No check here.) 
C
C     METHOD  C*1        I      (Uppercase) Type of fit to be used:
C                               'M' means Monotonic piecewise cubics;
C                               'B' means non-monotonic "Bessel"-type
C                                   piecewise cubics (looser fit);
C                               'L' means piecewise Linear fit
C
C     OFFSET             I      Value assigned to AREAS(1);
C                               added to all partial integrals
C
C     AREAS   (NDATA)      O    Integrals of Y spline up to each knot.
C                               AREAS(1) = OFFSET;
C                               AREAS(I) = OFFSET + integral from X(1) to X(I).
C
C  Significant local variables:
C
C     H, DEL         Delta X and forward difference derivative arrays.
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C  Procedures:
C
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C  Environment:  Fortran 90
C
C  History:
C
C     06/10/92  D.A.Saunders  LCSQUAD  adapted from LCSFIT.
C     06/12/02   "    "    "  LCSAREAS adapted from LCSQUAD and CSQUAD.
C                             Take advantage of marching through the
C                             data intervals by reusing the various
C                             2- and 3-point slope calculations,
C                             at the expense of repeated code.
C
C  Author:  David Saunders, ELORET/NASA Ames, Moffett Field, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA)
      CHARACTER, INTENT (IN)  :: METHOD * 1
      REAL,      INTENT (IN)  :: OFFSET
      REAL,      INTENT (OUT) :: AREAS (NDATA)

      REAL, PARAMETER ::
     &   ONE    = 1.0E+0,
     &   TWO    = 2.0E+0,
     &   THREE  = 3.0E+0,
     &   HALF   = ONE / TWO,
     &   THIRD  = ONE / THREE,
     &   FOURTH = ONE / (TWO + TWO)

C     Local variables:

      LOGICAL
     &   MONO
      INTEGER
     &   L
      REAL
     &   B (0:1), C, DELY (-1:1), D, DX, H (-1:1)

C     Procedures:

      REAL, EXTERNAL :: BESSEL, BRODLIE, BUTLAND, THREEPT

C     Execution:
C     ----------

C     Accumulate integrals for the cubics in successive intervals.

      AREAS (1) = OFFSET

      IF (METHOD .EQ. 'L' .OR. NDATA == 2) THEN ! Trapezoidal rule

         DO L = 2, NDATA
            AREAS (L) = (Y (L) + Y (L-1)) * HALF *
     &                  (X (L) - X (L-1)) + AREAS (L-1)
         END DO

      ELSE ! NDATA >= 3

         MONO = METHOD .EQ. 'M'

C        Left-most interval:

         H (0)    =  X (2) - X (1)
         DELY (0) = (Y (2) - Y (1)) / H (0) ! 2-pt. slopes
         H (1)    =  X (3) - X (2)
         DELY (1) = (Y (3) - Y (2)) / H (1)

         IF (MONO) THEN                     ! 3-pt. slopes
            B (0) = BUTLAND (0, H, DELY)
            B (1) = BRODLIE (1, H, DELY)
         ELSE
            B (0) = THREEPT (0, H, DELY)
            B (1) = BESSEL  (1, H, DELY)
         END IF

         DX = H (0)
         C  = (THREE * DELY (0) - TWO * B (0) - B (1)) / DX
         D  = ( -TWO * DELY (0) + B (0) + B (1)) / DX ** 2
         AREAS (2) = DX * (Y (1) + DX * (HALF * B (0) +
     &               DX * (THIRD * C + DX * FOURTH * D))) + AREAS (1)

C        Interior intervals:

         DO L = 3, NDATA - 1

            H (-1)    =  H (0)  ! Avoid recalculating earlier differences
            H ( 0)    =  H (1)
            H ( 1)    =  X (L+1) - X (L)
            DELY (-1) = DELY (0)
            DELY ( 0) = DELY (1)
            DELY ( 1) = (Y (L+1) - Y (L)) / H (1)

            IF (MONO) THEN
               B (0) = BRODLIE (0, H, DELY)
               B (1) = BRODLIE (1, H, DELY)
            ELSE
               B (0) = BESSEL  (0, H, DELY)
               B (1) = BESSEL  (1, H, DELY)
            END IF

            DX = H (0)
            C  = (THREE * DELY (0) - TWO * B (0) - B (1)) / DX
            D  = ( -TWO * DELY (0) + B (0) + B (1)) / DX ** 2
            AREAS (L) = DX * (Y (L-1) + DX * (HALF * B (0) +
     &                  DX * (THIRD * C + DX * FOURTH * D))) +
     &                  AREAS (L-1)
         END DO

C        Final interval:

         H (-1)    =  H (0)
         H ( 0)    =  H (1)
         DELY (-1) = DELY (0)
         DELY ( 0) = DELY (1)

         IF (MONO) THEN
            B (0) = BRODLIE (0, H, DELY)
            B (1) = BUTLAND (1, H, DELY)
         ELSE
            B (0) = BESSEL  (0, H, DELY)
            B (1) = THREEPT (1, H, DELY)
         END IF

         DX = H (0)
         C  = (THREE * DELY (0) - TWO * B (0) - B (1)) / DX
         D  = ( -TWO * DELY (0) + B (0) + B (1)) / DX ** 2
         L  = NDATA
         AREAS (L) = DX * (Y (L-1) + DX * (HALF * B (0) +
     &               DX * (THIRD * C + DX * FOURTH * D))) + AREAS (L-1)

      END IF

      END SUBROUTINE LCSAREAS
