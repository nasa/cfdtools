C+----------------------------------------------------------------------
C
      FUNCTION BUTLAND (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula, adjusted
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a modified forward or backward
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives, and the differencing direction is controlled by a flag.
C     See PLSFIT for more details, or THREEPT for the pure 3-pt. formula.
C
C        The "shape preserving adjustments" are from PCHIP, a monotone
C     piecewise cubic interpolation package by F. N. Fritsch.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BUTLAND R                 O    The function value is the adjusted
C                                    derivative.
C
C     Environment:  Fortran 90
C     ------------
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding, as THREEPT.
C     20 June 1991    DAS    Monotonic form renamed BUTLAND; THREEPT
C                            is now the pure 3-point formula.
C     04 Dec. 2002     "     SIGN work-arounds for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BUTLAND

C     Local constants.

      REAL, PARAMETER ::
     &   ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   THREE = 3.0E+0

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   DMAX, WEIGHT
      LOGICAL
     &   CONSTRAIN

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

C***  Avoid zero as the second argument of SIGN: Intel compiler misbehaves.

      IF (DEL (0) == ZERO) THEN

         BUTLAND = ZERO

      ELSE ! BUTLAND below cannot be zero

         WEIGHT  = -H (0) / (H (0) + H (STEP))
         BUTLAND = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C        Shape-preserving adjustments.  Note that we try to avoid overflow
C        by not multiplying quantities directly.

         IF (SIGN (ONE, BUTLAND) /= SIGN (ONE, DEL (0))) THEN

C           Defer to the estimate closest to the boundary.

            BUTLAND = ZERO

C******  ELSE IF (SIGN (ONE, DEL (0)) .NE. SIGN (ONE, DEL (STEP))) THEN
         ELSE

            IF (DEL (STEP) == ZERO) THEN
               CONSTRAIN = DEL (0) < ZERO
            ELSE
               CONSTRAIN = SIGN (ONE, DEL (0)) /= SIGN (ONE, DEL (STEP))
            END IF

            IF (CONSTRAIN) THEN

C              If the monotonicity switches, may need to bound the estimate.

               DMAX = THREE * DEL (0)
               IF (ABS (BUTLAND) > ABS (DMAX)) BUTLAND = DMAX
            END IF

         END IF

      END IF

C     Termination.
C     ------------

      END
