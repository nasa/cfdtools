C+----------------------------------------------------------------------
C
      SUBROUTINE LCSQUAD (NDATA, X, Y, XA, XB, METHOD, AREA)
C
C  One-liner:  Storage-efficient integration of discrete data
C
C  Description:
C
C        LCSQUAD estimates the integral on [XA, XB] defined by a discrete
C     distribution, using the storage-efficient spline interpolation of
C     LCSFIT.
C
C        Apart from giving results which are smoother with respect to
C     data perturbations than (say) the trapezoidal rule, the monotonic
C     option available for the underlying curve fit guards against possibly
C     wild excursions which may affect other polynomial-based methods.
C        
C  Arguments:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct and monotonic,
C                                    either ascending or descending.
C                                    (No check here.) 
C
C     XA,     R               I      Limits of integration, in either
C     XB                             order, possibly different from the
C                                    data abscissas.  The curve fit will
C                                    extrapolate if necessary, so beware.
C
C     METHOD  C*1             I      (Uppercase) Type of fit to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means Cyclic (periodic) end
C                                        conditions: loose fit assumed.
C
C     AREA    R                 O    Estimate of the integral of Y between
C                                    XA and XB. 
C
C  Significant local variables:
C
C     H, DEL         Delta X and forward difference derivative arrays.
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C  Procedures:
C
C     INTERVAL  1-D "interpolation" search utility.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C  Environment:  Fortran 90
C
C  History:
C
C     06/06/92  D.A.Saunders  Initial version using LCSFIT & QUANC8RC.
C     06/10/92    "    "      Adapted LCSFIT to do it more directly.
C     04/03/97    "    "      J. Reuther found pitched wing sections
C                             gave the wrong sign for ARROW. Now use
C                             X(NDATA) - X(1) instead of X(2) - X(1)
C                             to help wing fuel volume calculations.
C     05/06/98    "    "      Minor Fortran 90 updates.
C
C     Author: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA), XA, XB
      CHARACTER, INTENT (IN)  :: METHOD * 1
      REAL,      INTENT (OUT) :: AREA

      REAL, PARAMETER ::
     &   ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   TWO    = 2.0E+0,
     &   THREE  = 3.0E+0,
     &   HALF   = ONE / TWO,
     &   THIRD  = ONE / THREE,
     &   FOURTH = ONE / (TWO + TWO)

C     Local variables:

      LOGICAL
     &   CYCLIC, LINEAR, MONO
      INTEGER
     &   J, K, LEFT, LEFTA, LEFTB, RIGHT
      REAL
     &   ARROW, B (0:1), C, DELY (-1:1), D, H (-1:1), SUM, XLEFT, XR,
     &   XRIGHT

C     Procedures:

      REAL, EXTERNAL :: BESSEL, BRODLIE, BUTLAND, THREEPT

C     Statement function:

C     (This is the integral of the cubic for the interval
C     [X (LEFT), X (LEFT+1)] from X (LEFT) to X (LEFT) + DX.)

      REAL
     &   QUAD, DX
      QUAD (DX) = DX * (Y (LEFT) + DX * (HALF * B (0) +
     &            DX * (THIRD * C + DX * FOURTH * D)))

C     Execution:
C     ----------

      MONO   = METHOD .EQ. 'M'
      CYCLIC = METHOD .EQ. 'C'
      LINEAR = METHOD .EQ. 'L'

      IF (CYCLIC) THEN
         IF (Y (NDATA) /= Y (1)) STOP 'LCSQUAD: End pts. must match.'
      END IF

      ARROW = SIGN (ONE, X (NDATA) - X (1))

C     Locate the intervals containing XA and XB:

      LEFTA = 1
      CALL INTERVAL (NDATA, X, XA, ARROW, LEFTA)

      LEFTB = NDATA - 1
      CALL INTERVAL (NDATA, X, XB, ARROW, LEFTB)

C     Ensure that we traverse the X intervals in ascending order
C     (else we get the partial-interval handling messed up):

      IF (LEFTA <= LEFTB) THEN ! Reuse ARROW for sign of integral
         ARROW  = ONE
         XLEFT  = XA
         XRIGHT = XB
      ELSE
         ARROW  = -ONE
         LEFT   = LEFTA
         LEFTA  = LEFTB
         LEFTB  = LEFT
         XLEFT  = XB
         XRIGHT = XA
      END IF

C     Accumulate integrals for the cubics in successive intervals.
C     (Most of this loop is lifted from LCSFIT.)

      SUM = ZERO

      DO LEFT = LEFTA, LEFTB

         RIGHT = LEFT + 1

C        Compute derivatives by finite-differences:
C        ------------------------------------------

         IF (NDATA > 2 .AND. .NOT. LINEAR) THEN

            IF (LEFT == 1) THEN

               H (0) = X (2) - X (1)
               DELY (0) = (Y (2) - Y (1)) / H (0)
               H (1) = X (3) - X (2)
               DELY (1) = (Y (3) - Y (2)) / H (1)

               IF (CYCLIC) THEN ! Loose fit assumed
                  H (-1) = X (NDATA) - X (NDATA - 1)
                  DELY (-1) = (Y (NDATA) - Y (NDATA - 1)) / H (-1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE
                  IF (MONO) THEN
                     B (0) = BUTLAND (0, H, DELY)
                     B (1) = BRODLIE (1, H, DELY)
                  ELSE
                     B (0) = THREEPT (0, H, DELY)
                     B (1) = BESSEL  (1, H, DELY)
                  END IF
               END IF

            ELSE IF (RIGHT == NDATA) THEN

               H(-1) = X (LEFT) - X (LEFT-1)
               DELY(-1) = (Y (LEFT) - Y (LEFT-1)) / H (-1)
               H (0) = X (RIGHT) - X (LEFT)
               DELY (0) = (Y (RIGHT) - Y (LEFT))  / H (0)

               IF (CYCLIC) THEN
                  H (1) = X (2) - X (1)
                  DELY (1) = (Y (2) - Y (1)) / H (1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE

                  IF (MONO) THEN
                     B (0) = BRODLIE (0, H, DELY)
                     B (1) = BUTLAND (1, H, DELY)
                  ELSE
                     B (0) = BESSEL  (0, H, DELY)
                     B (1) = THREEPT (1, H, DELY)
                  END IF
               END IF

            ELSE

               K = LEFT
               DO J = -1, +1
                  H (J)    =  X (K) - X (K-1)
                  DELY (J) = (Y (K) - Y (K-1)) / H (J)
                  K = K + 1
               END DO

C              Select interpolation scheme:
C              ----------------------------

C              Compute (possibly adjusted) first derivatives at both
C              left- and right-hand endpoints of the interval:

               IF (MONO) THEN

C                 Monotone - use Brodlie modification of Butland's
C                 formula to adjust the derivatives at the knots.

                  B (0) = BRODLIE (0, H, DELY)
                  B (1) = BRODLIE (1, H, DELY)

               ELSE ! METHOD = 'B'

C                 Bessel - use central difference formula at the knots.

                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)

               END IF

            END IF

C           Compute the remaining cubic coefficients.

            C = (THREE * DELY (0) - TWO * B (0) - B (1)) / H (0)
            D = ( -TWO * DELY (0) + B (0) + B (1)) / H (0) ** 2

         ELSE ! NDATA = 2 .OR. METHOD = 'L'

C           Degenerate case (linear):
C           -------------------------

            B (0) = (Y (RIGHT) - Y (LEFT)) / (X (RIGHT) - X (LEFT))
            C     = ZERO
            D     = ZERO

         END IF

C        Evaluate the integral of the cubic for this interval:
C        -----------------------------------------------------

C        (Make end adjustments in case XA or XB is not at an X (I).)

         IF (LEFT == LEFTA) SUM = SUM - QUAD (XLEFT - X (LEFT))

         IF (LEFT /= LEFTB) THEN
            XR = X (RIGHT)  ! For complete interior intervals
         ELSE
            XR = XRIGHT     ! No matter where XB is
         END IF

         SUM = SUM + QUAD (XR - X (LEFT))

      END DO
   
      AREA = SUM * ARROW

      RETURN

      END SUBROUTINE LCSQUAD
