C+----------------------------------------------------------------------
C
      FUNCTION BRODLIE (J, H, DEL)
C
C     One-liner: First derivative, adjusted for monotonicity
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        BRODLIE is intended to be used by PLSFIT for determining end
C     conditions on an interval for monotonic interpolation by piecewise
C     cubics. The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives. See the PLSFIT header for more details.
C
C        The method is due to Brodlie, Butland, Carlson, and Fritsch,
C     as referenced in the PLSFIT header.
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
C     BRODLIE R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     04 Dec. 2002    DAS    SIGN work-around for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE, THIRD
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   THIRD  = ONE / 3.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   BRODLIE, H (-1:1), DEL (-1:1)

C     Local variables.

      REAL
     &   ALPHA, PRODUCT

C     Execution.
C     ----------

C     Compare the algebraic signs of the two DEL's.  Have to test that
C     at least one is positive to avoid a zero denominator (this fancy
C     test permits one term to be zero, but the answer below is zero
C     anyway in these cases).  The trick is to work around the SIGN
C     function, which returns positive even if its 2nd argument is zero.

C**** NO:  SIGN misbehaves on Intel systems when the 2nd argument is zero.

      PRODUCT = DEL (J - 1) * DEL (J)

      IF (PRODUCT == ZERO) THEN

         BRODLIE = ZERO

      ELSE IF (SIGN (ONE, -DEL (J - 1)) .NE. SIGN (ONE, DEL (J))) THEN

C        Form "weighted harmonic mean" of the two finite-difference
C        derivative approximations.  Note that we try to avoid overflow
C        by not multiplying them together directly.

         ALPHA   = THIRD * (ONE + H (J) / (H (J - 1) + H (J)))
         BRODLIE = PRODUCT / (ALPHA * DEL (J) +
     &      (ONE - ALPHA) * DEL (J - 1))
      ELSE

C        The signs differ, so make this point a local extremum.

         BRODLIE = ZERO
      END IF

C     Termination.
C     ------------

      RETURN
      END
